program my_rte_rrtmgp_allsky
  use, intrinsic :: iso_fortran_env, & 
                             only: output_unit
  use mo_rte_kind,           only: wp, i8, wl
  use mo_optical_props,      only: ty_optical_props, &
                                   ty_optical_props_arry, ty_optical_props_1scl, ty_optical_props_2str
  use mo_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp
  use mo_cloud_optics_rrtmgp,only: ty_cloud_optics_rrtmgp
  use mo_aerosol_optics_rrtmgp_merra ! Includes aerosol type integers
  use mo_gas_concentrations, only: ty_gas_concs
  use mo_source_functions,   only: ty_source_func_lw
  use mo_fluxes,             only: ty_fluxes_broadband
  use mo_rte_lw,             only: rte_lw
  use mo_rte_sw,             only: rte_sw
  use mo_load_coefficients,  only: load_and_init
  use mo_load_cloud_coefficients, &
                             only: load_cld_lutcoeff, load_cld_padecoeff
  use mo_load_aerosol_coefficients, &
                             only: load_aero_lutcoeff
  use mo_rte_config,         only: rte_config_checks
  use mo_load_profiles,      only: load_profiles, load_vert_grid_info, load_2d_sw, load_2d_lw
  use mo_gas_optics_constants,   only: grav
  use mo_gas_optics_util_string, only: string_loc_in_array
  implicit none
  ! ----------------------------------------------------------------------------------
  ! Variables
  ! ----------------------------------------------------------------------------------
  ! Arrays: dimensions (col, lay)
  real(wp), dimension(:,:),   allocatable :: p_lay, t_lay, p_lev, t_lev, h_lev, dh_lev ! t_lev is only needed for LW
  real(wp), dimension(:,:),   allocatable :: q, o3, col_dry, qc, qi, n2o, ch4
  real(wp), dimension(:,:),   allocatable :: q_mmr, o3_mmr, qc_mmr, qi_mmr, n2o_mmr, ch4_mmr
  real(wp), dimension(:,:),   allocatable :: temp_array

  !
  ! Longwave only
  !
  real(wp), dimension(:),     allocatable :: t_sfc, solin
  real(wp), dimension(:,:),   allocatable :: emis_sfc ! First dimension is band
  !
  ! Shortwave only
  !
  real(wp), dimension(:),     allocatable :: mu0
  real(wp), dimension(:,:),   allocatable :: sfc_alb_dir, sfc_alb_dif ! First dimension is band
  !
  ! Gas concentrations 
  !
  character(len=3), dimension(8), parameter :: &
                   gas_names = ['h2o', 'co2', 'o3 ', 'n2o', 'co ', 'ch4', 'o2 ', 'n2 ']
  !
  ! Source functions
  !
  !   Longwave
  type(ty_source_func_lw), save               :: lw_sources
  !   Shortwave
  real(wp), dimension(:,:), allocatable, save :: toa_flux
  !
  ! Clouds
  !
  real(wp), allocatable, dimension(:,:) :: lwp, iwp, rel, rei
  logical,  allocatable, dimension(:,:) :: cloud_mask
  !
  ! Aerosols
  !
  logical :: cell_has_aerosols
  integer,  dimension(:,:), allocatable :: aero_type 
                                           ! MERRA2/GOCART aerosol type
  real(wp), dimension(:,:), allocatable :: aero_size
                                           ! Aerosol size for dust and sea salt
  real(wp), dimension(:,:), allocatable :: aero_mass
                                           ! Aerosol mass column (kg/m2)
  real(wp), dimension(:,:), allocatable :: relhum
                                           ! Relative humidity (fraction)
  logical, dimension(:,:), allocatable  :: aero_mask
                                           ! Aerosol mask

  !
  ! Output variables
  !
  real(wp), dimension(:,:), target, &
                            allocatable :: flux_up, flux_dn, flux_dir
  !
  ! Derived types from the RTE and RRTMGP libraries
  !
  type(ty_gas_optics_rrtmgp)   :: k_dist
  type(ty_cloud_optics_rrtmgp) :: cloud_optics
  type(ty_aerosol_optics_rrtmgp_merra)   & 
                               :: aerosol_optics
  type(ty_gas_concs)           :: gas_concs, gas_concs_garand, gas_concs_1col
  class(ty_optical_props_arry), &
                 allocatable   :: atmos, clouds, aerosols
  type(ty_fluxes_broadband)    :: fluxes

  !
  ! Inputs to RRTMGP
  !
  logical :: top_at_1, is_sw, is_lw

  integer  :: nbnd, ngpt
  integer  :: icol, ilay, ibnd, iloop, igas

  character(len=8) :: char_input
  integer :: nUserArgs, nloops, ncol, nlay
  ! logical :: write_fluxes = .false.
  logical :: do_clouds = .false., use_luts = .true. 
  logical :: do_aerosols = .false.
  integer, parameter :: ngas = 8

  character(len=256) :: output_file, k_dist_file, cloud_optics_file, aerosol_optics_file
  !
  ! Timing variables
  !
  integer(kind=i8)              :: start, finish, start_all, finish_all, clock_rate
  real(wp)                      :: avg, mint
  integer(kind=i8), allocatable :: elapsed(:)
  ! my vars
  character(len=256) :: cs_mli_file, cs_grid_info_file
  real(wp), dimension(:),   allocatable :: hyai, hybi, hyam, hybm, lwup, ps
  real(wp), allocatable :: p0
  real(wp), parameter :: sigma = 5.670374419e-8_wp ! Stefan-Boltzmann constant 
  !----------
  ! NAR OpenMP CPU directives in compatible with OpenMP GPU directives
  !!$omp threadprivate( lw_sources, toa_flux, flux_up, flux_dn, flux_dir )
  ! ----------------------------------------------------------------------------------
  ! Code
  ! ----------------------------------------------------------------------------------
  !
  ! Parse command line: rrtmgp_allsky ncol nlay nreps kdist [clouds [aerosols]] 

  !
  nUserArgs = command_argument_count()
  if (nUserArgs <  7) call stop_on_err("Usage: rrtmgp_allsky ncol nlay nreps output_file gas-optics climsim_lowres_grid_file climsim_mli_file [cloud-optics [aerosol-optics]]")

  call get_command_argument(1, char_input)
  read(char_input, '(i8)') ncol
  if(ncol <= 0) call stop_on_err("Specify positive ncol.")

  call get_command_argument(2, char_input)
  read(char_input, '(i8)') nlay
  if(nlay <= 0) call stop_on_err("Specify positive nlay.")

  call get_command_argument(3, char_input)
  read(char_input, '(i8)') nloops
  if(nloops <= 0) call stop_on_err("Specify positive nreps (number of times to do ncol examples.")

  call get_command_argument(4,output_file)
  call get_command_argument(5,k_dist_file)

  call get_command_argument(6,cs_grid_info_file)
  call get_command_argument(7,cs_mli_file)

  if (nUserArgs >= 8) then 
    call get_command_argument(8,cloud_optics_file)
    do_clouds = .true. 
  end if 
  if (nUserArgs >= 9) then 
    call get_command_argument(9,aerosol_optics_file)
    do_aerosols = .true. 
  end if 
  if (nUserArgs >  9) print *, "Ignoring command line arguments beyond the first seven..."
  ! -----------------------------------------------------------------------------------
  allocate(p_lay(ncol, nlay), t_lay(ncol, nlay), p_lev(ncol, nlay+1), t_lev(ncol, nlay+1), h_lev(ncol, nlay+1), dh_lev(ncol, nlay))
  allocate(q    (ncol, nlay),    o3(ncol, nlay), qc   (ncol, nlay),   qi   (ncol, nlay), n2o(ncol, nlay), ch4(ncol, nlay), ps(ncol))
  allocate(q_mmr(ncol, nlay),o3_mmr(ncol, nlay), qc_mmr(ncol, nlay), qi_mmr(ncol, nlay), n2o_mmr(ncol, nlay), ch4_mmr(ncol, nlay))
  !$acc        data create(   p_lay, t_lay, p_lev, t_lev, q, o3)
  !$omp target data map(alloc:p_lay, t_lay, p_lev, t_lev, q, o3)
  ! ------------------- my stuff ------------------
  ! cs_grid_info_file = "/p/scratch/icon-a-ml/heuer1/LEAP/ClimSim_low-res/ClimSim_low-res_grid-info.nc"
  allocate(hyai(nlay+1), hybi(nlay+1), hyam(nlay), hybm(nlay), p0)
  call load_vert_grid_info(cs_grid_info_file, ncol, nlay, hyai, hybi, hyam, hybm, p0)
  ! print *, p0
  ! print *, hyai
  ! print *, hybi
  ! print *, hyam
  ! print *, hybm
  ! print *, cs_grid_info_file
  ! print *, cs_mli_file

  ! cs_mli_file = "/p/scratch/icon-a-ml/heuer1/LEAP/ClimSim_low-res/train/0001-02/E3SM-MMF.mli.0001-02-01-00000.nc"
  call load_profiles(cs_mli_file,&
                     ncol,&
                     nlay,&
                     nbnd,&
                     hyai,&
                     hybi,&
                     hyam,&
                     hybm,&
                     p0,&
                     ps,&
                     p_lay,&
                     t_lay,&
                     p_lev,&
                     q_mmr,&
                     qc_mmr,&
                     qi_mmr,&
                     o3_mmr,&
                     n2o_mmr,&
                     ch4_mmr) 
        
  ! from https://www.weather.gov/media/epz/wxcalc/pressureAltitude.pdf
  do ilay=1,nlay+1
    do icol=1,ncol
      h_lev(icol,ilay) = 0.3048_wp*145366.45_wp*(1 - (p_lev(icol,ilay)/ps(icol))**0.190284)
    end do
  end do
  do ilay=1,nlay
    do icol=1,ncol
      dh_lev(icol,ilay) = h_lev(icol,ilay) - h_lev(icol,ilay+1)
    end do
  end do
  ! print *, size(h_lev, 1)
  ! print *, "------------"
  ! print *, sum(h_lev, 1)/size(h_lev, 1)
  ! print *, "------------"
  ! print *, sum(dh_lev, 1)/size(dh_lev, 1)
  ! print *, p_lay(192,:)
  ! print *, "------------------"
  ! print *, p_lev(192,:)
  ! print *, "------------------"
  ! print *, t_lay(192,:)
  ! print *, "------------------"
  
  call get_gas_vmr(ncol, nlay, q_mmr, 'H2O', q)
  call get_gas_vmr(ncol, nlay, qc_mmr, 'H2O', qc)
  call get_gas_vmr(ncol, nlay, qi_mmr, 'H2O', qi)
  call get_gas_vmr(ncol, nlay, o3_mmr, 'O3', o3)
  call get_gas_vmr(ncol, nlay, n2o_mmr, 'N2O', n2o)
  call get_gas_vmr(ncol, nlay, ch4_mmr, 'CH4', ch4)
  ! ------------------- my stuff ------------------
  ! call compute_profiles(300._wp, ncol, nlay, p_lay, t_lay, p_lev, t_lev, q, o3)
  ! print *, p_lev(1,:)
  ! print *, p_lay(1,:)
  ! print *, (p_lev(1,1:60) + p_lev(1,2:61))*0.5

  call stop_on_err(gas_concs%init(gas_names))
  call stop_on_err(gas_concs%set_vmr("h2o", q )) 
  call stop_on_err(gas_concs%set_vmr("o3",  o3)) 
  call stop_on_err(gas_concs%set_vmr("co2", 388.717e-6_wp)) ! from &chem_surfvals_nl in atm_in (e3sm)
  call stop_on_err(gas_concs%set_vmr("ch4", ch4)) 
  call stop_on_err(gas_concs%set_vmr("n2o", n2o)) 
  ! call stop_on_err(gas_concs%set_vmr("co2", 348.e-6_wp)) 
  ! call stop_on_err(gas_concs%set_vmr("ch4", 1650.e-9_wp)) 
  ! call stop_on_err(gas_concs%set_vmr("n2o", 306.e-9_wp)) 
  call stop_on_err(gas_concs%set_vmr("n2",  0.7808_wp)) 
  call stop_on_err(gas_concs%set_vmr("o2",  0.2095_wp)) 
  call stop_on_err(gas_concs%set_vmr("co",  0._wp)) 
  ! ----------------------------------------------------------------------------
  ! load data into classes
  call load_and_init(k_dist, k_dist_file, gas_concs)
  is_sw = k_dist%source_is_external()
  is_lw = .not. is_sw
  if (do_clouds) then 
    !
    ! Should also try with Pade calculations
    !  call load_cld_padecoeff(cloud_optics, cloud_optics_file)
    !
    if(use_luts) then
      call load_cld_lutcoeff (cloud_optics, cloud_optics_file)
    else
      call load_cld_padecoeff(cloud_optics, cloud_optics_file)
    end if
    call stop_on_err(cloud_optics%set_ice_roughness(2))
  end if

  if (do_aerosols) then 
    !
    ! Load aerosol optics coefficients from lookup tables
    !
    call load_aero_lutcoeff (aerosol_optics, aerosol_optics_file)
  end if 

  ! ----------------------------------------------------------------------------
  !
  ! Problem sizes
  !
  nbnd = k_dist%get_nband()
  ! print *,nbnd
  ngpt = k_dist%get_ngpt()
  top_at_1 = p_lay(1, 1) < p_lay(1, nlay)

  ! ----------------------------------------------------------------------------
  ! LW calculations neglect scattering; SW calculations use the 2-stream approximation
  !   Here we choose the right variant of optical_props.
  !
  if(is_sw) then
    allocate(ty_optical_props_2str::atmos)
  else
    allocate(ty_optical_props_1scl::atmos)
  end if
  !
  ! Allocate arrays for the optical properties themselves.
  !
  select type(atmos)
    class is (ty_optical_props_1scl)
      !$acc enter data copyin(atmos)
      call stop_on_err(atmos%alloc_1scl(ncol, nlay, k_dist))
      !$acc enter data copyin(atmos) create(atmos%tau)
      !$omp target enter data map(alloc:atmos%tau)
    class is (ty_optical_props_2str)
      call stop_on_err(atmos%alloc_2str( ncol, nlay, k_dist))
      !$acc enter data copyin(atmos) create(atmos%tau, atmos%ssa, atmos%g)
      !$omp target enter data map(alloc:atmos%tau, atmos%ssa, atmos%g)
    class default
      call stop_on_err("my_rte_rrtmgp_allsky: Don't recognize the kind of optical properties ")
  end select

  ! ----------------------------------------------------------------------------
  !  Boundary conditions depending on whether the k-distribution being supplied
  !   is LW or SW
  if(is_sw) then
    ! toa_flux is threadprivate
    !!$omp parallel
    allocate(toa_flux(ncol, ngpt))
    allocate(sfc_alb_dir(nbnd, ncol), sfc_alb_dif(nbnd, ncol), mu0(ncol), solin(ncol))
    call load_2d_sw(cs_mli_file,&
                      ncol,&
                      nlay,&
                      nbnd,&
                      k_dist%get_band_lims_wavenumber(),&
                      solin,&
                      mu0,&
                      sfc_alb_dir,&
                      sfc_alb_dif)

    ! print *, solin(192)
    ! print *, "------------------"
    ! print *, mu0(192)
    ! print *, "------------------"
    ! print *, sfc_alb_dir(:,192)
    ! print *, "------------------"
    ! print *, sfc_alb_dif(:,192)
    ! print *, "------------------"

    !!$omp end parallel
    !
    !$acc         enter data create(   sfc_alb_dir, sfc_alb_dif, mu0)
    !$omp target  enter data map(alloc:sfc_alb_dir, sfc_alb_dif, mu0)
    ! Ocean-ish values for no particular reason
    !$acc kernels
    !$omp target
    ! sfc_alb_dir = 0.06_wp
    ! sfc_alb_dif = 0.06_wp
    ! mu0 = .86_wp
    !$acc end kernels
    !$omp end target
  else
    ! lw_sorces is threadprivate
    !!$omp parallel
    call stop_on_err(lw_sources%alloc(ncol, nlay, k_dist))
    !!$omp end parallel

    allocate(t_sfc(ncol), emis_sfc(nbnd, ncol), lwup(ncol))
    call load_2d_lw(cs_mli_file,&
                      ncol,&
                      nlay,&
                      nbnd,&
                      lwup)
                    
    ! print *, lwup(192)
    ! print *, "------------------"

    !$acc         enter data create   (t_sfc, emis_sfc)
    !$omp target  enter data map(alloc:t_sfc, emis_sfc)
    ! Surface temperature
    !$acc kernels
    !$omp target
    ! emis_sfc = 0.98_wp
    ! Set to 1 instead after this note (to be discussed): https://github.com/E3SM-Project/ACME-ECP/blob/3912d6af813c80c6b7ca6264916cafc0a9e0069e/components/cam/src/physics/rrtmgp/radiation.F90#L1289-L1294
    emis_sfc = 1.0_wp

    ! t_sfc = t_lev(1, merge(nlay+1, 1, top_at_1))
    t_sfc = (lwup/(emis_sfc(1,:)*sigma))**(1.0/4.0)
    ! print *, sum(t_sfc)/size(t_sfc)

    !$acc end kernels
    !$omp end target
  end if
  ! ----------------------------------------------------------------------------
  !
  ! Fluxes
  !
  !!$omp parallel
  allocate(flux_up(ncol,nlay+1), flux_dn(ncol,nlay+1))
  !!$omp end parallel

  !$acc         data create(   flux_up, flux_dn)
  !$omp target  data map(alloc:flux_up, flux_dn)
  if(is_sw) then
    allocate(flux_dir(ncol,nlay+1))
    !$acc enter data create(flux_dir)
    !$omp target enter data map(alloc:flux_dir)
  end if

  if (do_clouds)   call compute_clouds
  if (do_aerosols) call compute_aerosols

  ! ----------------------------------------------------------------------------
  !
  ! Multiple iterations for big problem sizes, and to help identify data movement
  !   For CPUs we can introduce OpenMP threading over loop iterations
  !
  allocate(elapsed(nloops))
  !
  call system_clock(start_all)
  !
  !!$omp parallel do firstprivate(fluxes)
  do iloop = 1, nloops
    ! Omit the checks starting with the second iteration
    if (iloop > 1) call rte_config_checks(logical(.false., wl))

    call system_clock(start)
    !
    ! Cloud optics
    !
    if(do_clouds) & 
      call stop_on_err(cloud_optics%cloud_optics(lwp, iwp, rel, rei, clouds))
    !
    ! Aerosol optics
    !
    if(do_aerosols) & 
      call stop_on_err(aerosol_optics%aerosol_optics(aero_type, aero_size,  &
                                                     aero_mass, relhum, aerosols))
    !
    ! Solvers
    !
    fluxes%flux_up => flux_up(:,:)
    fluxes%flux_dn => flux_dn(:,:)
    if(is_lw) then
      !
      ! Should we allocate these once, rather than once per loop? They're big. 
      ! 
      !$acc        data create(   lw_sources, lw_sources%lay_source,     lw_sources%lev_source) &
      !$acc             create(               lw_sources%sfc_source,     lw_sources%sfc_source_Jac)
      !$omp target data map(alloc:            lw_sources%lay_source,     lw_sources%lev_source) &
      !$omp             map(alloc:            lw_sources%sfc_source,     lw_sources%sfc_source_Jac)
      call stop_on_err(k_dist%gas_optics(p_lay, p_lev, &
                                         t_lay, t_sfc, &
                                         gas_concs,    &
                                         atmos,        &
                                         lw_sources))!,   &
                                        !  tlev = t_lev))
      if(do_clouds)   call stop_on_err(clouds%increment(atmos))
      if(do_aerosols) call stop_on_err(aerosols%increment(atmos))
      call stop_on_err(rte_lw(atmos, top_at_1, &
                              lw_sources,      &
                              emis_sfc,        &
                              fluxes))
      !$acc        end data
      !$omp end target data

    else
      !$acc         data create(   toa_flux)
      !$omp target  data map(alloc:toa_flux)
      fluxes%flux_dn_dir => flux_dir(:,:)

      call stop_on_err(k_dist%gas_optics(p_lay, p_lev, &
                                         t_lay,        &
                                         gas_concs,    &
                                         atmos,        &
                                         toa_flux))
      ! rescale toa_flux
      do icol = 1, ncol 
        ! print *, solin(icol) / SUM(toa_flux(icol,:))
        toa_flux(icol,:) = toa_flux(icol,:) * solin(icol) / SUM(toa_flux(icol,:))
      end do
      ! end rescale toa_flux
      if(do_clouds) then 
        call stop_on_err(clouds%delta_scale())
        call stop_on_err(clouds%increment(atmos))
      end if 
      if(do_aerosols) then 
        call stop_on_err(aerosols%delta_scale())
        call stop_on_err(aerosols%increment(atmos))
      end if 
      call stop_on_err(rte_sw(atmos, top_at_1, &
                              mu0,   toa_flux, &
                              sfc_alb_dir, sfc_alb_dif, &
                              fluxes))
      ! print *, maxval(fluxes%flux_up(:,:), 1)
      ! print *, "------------------"
      ! print *, maxval(fluxes%flux_dn(:,:), 1)

      !$acc        end data   
      !$omp end target data
    end if
    call system_clock(finish, clock_rate)
    elapsed(iloop) = finish - start
  end do
  !
  call system_clock(finish_all, clock_rate)

  avg  = sum( elapsed(merge(2,1,nloops>1):) ) / real(merge(nloops-1,nloops,nloops>1))
  mint = minval(elapsed) 

  ! What to print? 
  !   ncol, nlay, ngpt; are clouds used, are aerosols used; time per column, total, min; 
  print *, " ncol   nlay   ngpt  clouds aerosols time_per_col_ms nloops time_total_s time_min_s"
  write(output_unit, '(3(i6, 1x), 6x, 2(i1, 8x), 1x, f7.3, 1x, i6, 2x, 2(4x,f7.3))') & 
    ncol, nlay, ngpt, merge(1,0,do_clouds), merge(1,0,do_aerosols),  & 
    avg/(real(ncol) * (1.0e-3*clock_rate)),  nloops,  sum(elapsed) / real(clock_rate),  mint / real(clock_rate)

  call write_fluxes

  ! 
  ! Memory for bounday conditions on the GPU was allocated with unstructured data dataments 
  !   (acc enter data). Deallocate it expliicity 
  !
  if(is_lw) then
    !$acc        exit data delete(     t_sfc, emis_sfc)
    !$omp target exit data map(release:t_sfc, emis_sfc)
  else
    !$acc        exit data delete(     sfc_alb_dir, sfc_alb_dif, mu0)
    !$omp target exit data map(release:sfc_alb_dir, sfc_alb_dif, mu0)
  end if
  
  !
  ! Clouds and aerosols also used enter data  
  !
  if(do_clouds) then
#ifndef _CRAYFTN
    ! ACCWA cloud_mask is already deallocated, Cray does not currently allow
    ! ACCWA delete in that case
    !$acc        exit data delete(     cloud_mask, lwp, iwp, rel, rei)
#endif
    !$omp target exit data map(release:cloud_mask, lwp, iwp, rel, rei)
    select type(clouds)
      class is (ty_optical_props_1scl)
        !$acc        exit data delete     (clouds%tau, clouds)
        !$omp target exit data map(release:clouds%tau)
      class is (ty_optical_props_2str)
        !$acc        exit data delete     (clouds%tau, clouds%ssa, clouds%g, clouds)
        !$omp target exit data map(release:clouds%tau, clouds%ssa, clouds%g)
    end select
    ! 
    ! Explicit finalization of cloud optical properties - not really necessary since memory 
    !   will be freed when the program ends, but useful for testing 
    !
    call clouds%finalize
  end if
  if(do_aerosols) then
    !$acc        exit data delete(     aero_type, aero_size, aero_mass, relhum)
    !$omp target exit data map(release:aero_type, aero_size, aero_mass, relhum)
    select type(aerosols)
      class is (ty_optical_props_1scl)
        !$acc        exit data delete     (aerosols%tau, aerosols)
        !$omp target exit data map(release:aerosols%tau)
      class is (ty_optical_props_2str)
        !$acc        exit data delete     (aerosols%tau, aerosols%ssa, aerosols%g, aerosols)
        !$omp target exit data map(release:aerosols%tau, aerosols%ssa, aerosols%g)
    end select
    ! 
    ! Explicit finalization of aerosol optical properties - not really necessary since memory 
    !   will be freed when the program ends, but useful for testing 
    !
    call aerosols%finalize
  end if
  !
  ! k-distribution
  !
  call k_dist%finalize
  
  if(.not. is_lw) then
    !$acc        exit data delete(     flux_dir)
    !$omp target exit data map(release:flux_dir)
  end if

  ! fluxes - but not flux_dir, which used enter data 
  !$acc end        data 
  !$omp end target data
  ! p_lay etc
  !$acc end        data 
  !$omp end target data
contains
  ! ----------------------------------------------------------------------------------
  subroutine stop_on_err(error_msg)
    use iso_fortran_env, only : error_unit
    character(len=*), intent(in) :: error_msg

    if(error_msg /= "") then
      write (error_unit,*) trim(error_msg)
      write (error_unit,*) "rrtmgp_allsky stopping"
      error stop 1
    end if
  end subroutine stop_on_err
  ! --------------------------------------------------------------------------------------
  !
  subroutine compute_clouds 
    real(wp) :: rel_val, rei_val
    ! for icon adaption
    real(wp) :: reimin, reimax, relmin, relmax, re_cryst_scal, re_drop_scal, ziwc, zlwc, zkap, effective_radius
    real(wp) :: zprat, zn1, zn2
    ! declare zcdnc as array for debugging, computation of zcdnc has to be shifted before if condition later
    ! real(wp), dimension(ncol,nlay) :: zcdnc
    real(wp) :: zcdnc

    real(wp), parameter :: cn1lnd = 20._wp
    real(wp), parameter :: cn2lnd = 180._wp
    real(wp), parameter :: cn1sea = 20._wp
    real(wp), parameter :: cn2sea = 80._wp
    ! mean cdnc profile analyzed from simulation data (changed to computation as in icon)
    ! real(wp) :: mean_cdnc(60) = (/2.00000000e+07, 2.00000000e+07, 2.00000000e+07, 2.00000000e+07,&
    !                               2.00000000e+07, 2.00000000e+07, 2.00000000e+07, 2.00000000e+07,&
    !                               2.00000000e+07, 2.00000000e+07, 2.00000000e+07, 2.00000000e+07,&
    !                               2.00000000e+07, 2.00000000e+07, 2.00000000e+07, 2.00000000e+07,&
    !                               2.00000000e+07, 2.00000000e+07, 2.00000000e+07, 2.00000000e+07,&
    !                               2.00000000e+07, 2.00000000e+07, 2.00000000e+07, 2.00000000e+07,&
    !                               2.00000000e+07, 2.00000000e+07, 2.00000000e+07, 2.00000000e+07,&
    !                               2.00000000e+07, 2.00000000e+07, 2.00000000e+07, 2.00000000e+07,&
    !                               2.00000000e+07, 2.00000000e+07, 2.00000000e+07, 2.00000680e+07,&
    !                               2.00023200e+07, 2.00295560e+07, 2.01922880e+07, 2.07834940e+07,&
    !                               2.22821960e+07, 2.51960240e+07, 2.98459080e+07, 3.62477360e+07,&
    !                               4.41368320e+07, 5.30814960e+07, 6.26067720e+07, 7.22813280e+07,&
    !                               8.17591280e+07, 9.07884320e+07, 9.92019200e+07, 1.04135896e+08,&
    !                               1.05483576e+08, 1.06020616e+08, 1.06287360e+08, 1.06438648e+08,&
    !                               1.06539488e+08, 1.06609136e+08, 1.06655440e+08, 1.06689352e+08/)
      ! adapted from ICON radiation code
      real(wp), parameter :: droplet_scale = 1.0e2
      real(wp), parameter :: pi = acos(-1._wp)
      real(wp), parameter :: rhoh2o = 1000._wp
      real(wp), parameter :: ccwmin = 1.e-7_wp      ! min condensate for lw cloud opacity
      real(wp), parameter :: zkap_cont = 1.143_wp   ! continental (Martin et al. ) breadth param
      real(wp), parameter :: zkap_mrtm = 1.077_wp    ! maritime (Martin et al.) breadth parameter
    ! 
    ! Variable and memory allocation 
    !
    if(is_sw) then
      allocate(ty_optical_props_2str::clouds)
    else
      allocate(ty_optical_props_1scl::clouds)
    end if
    ! Clouds optical props are defined by band
    call stop_on_err(clouds%init(k_dist%get_band_lims_wavenumber()))

    select type(clouds)
      class is (ty_optical_props_1scl)
        call stop_on_err(clouds%alloc_1scl(ncol, nlay))
        !$acc enter data copyin(clouds) create(clouds%tau)
        !$omp target enter data map(alloc:clouds%tau)
      class is (ty_optical_props_2str)
        call stop_on_err(clouds%alloc_2str(ncol, nlay))
        !$acc enter data copyin(clouds) create(clouds%tau, clouds%ssa, clouds%g)
        !$omp target enter data map(alloc:clouds%tau, clouds%ssa, clouds%g)
      class default
        call stop_on_err("my_rte_rrtmgp_allsky: Don't recognize the kind of optical properties ")
    end select
    !
    ! Cloud physical properties 
    !
    allocate(lwp(ncol,nlay), iwp(ncol,nlay), &
             rel(ncol,nlay), rei(ncol,nlay), cloud_mask(ncol,nlay))
    !$acc enter        data create(   cloud_mask, lwp, iwp, rel, rei)
    !$omp target enter data map(alloc:cloud_mask, lwp, iwp, rel, rei)

    ! ----------------------------------- default method ------------------------------------ !
    ! ! Restrict clouds to troposphere (> 100 hPa = 100*100 Pa)
    ! !   and not very close to the ground (< 900 hPa), and
    ! !   put them in 2/3 of the columns since that's roughly the
    ! !   total cloudiness of earth
    ! rel_val = 0.5 * (cloud_optics%get_min_radius_liq() + cloud_optics%get_max_radius_liq())
    ! rei_val = 0.5 * (cloud_optics%get_min_radius_ice() + cloud_optics%get_max_radius_ice())
    ! !$acc                         parallel loop    collapse(2) copyin(t_lay) copyout( lwp, iwp, rel, rei)
    ! !$omp target teams distribute parallel do simd collapse(2) map(to:t_lay) map(from:lwp, iwp, rel, rei)
    ! do ilay=1,nlay
    !   do icol=1,ncol
    !     cloud_mask(icol,ilay) = p_lay(icol,ilay) > 100._wp * 100._wp .and. &
    !                             p_lay(icol,ilay) < 900._wp * 100._wp .and. &
    !                             mod(icol, 3) /= 0
    !     !
    !     ! Ice and liquid will overlap in a few layers
    !     !
    !     ! lwp(icol,ilay) = merge(10._wp,  0._wp, cloud_mask(icol,ilay) .and. t_lay(icol,ilay) > 263._wp)
    !     ! iwp(icol,ilay) = merge(10._wp,  0._wp, cloud_mask(icol,ilay) .and. t_lay(icol,ilay) < 273._wp)

    !     lwp(icol,ilay) = (p_lev(icol,ilay+1) - p_lev(icol,ilay)) * qc(icol,ilay) / grav * 1e3_wp
    !     iwp(icol,ilay) = (p_lev(icol,ilay+1) - p_lev(icol,ilay)) * qi(icol,ilay) / grav * 1e3_wp

    !     rel(icol,ilay) = merge(rel_val, 0._wp, lwp(icol,ilay) > 0._wp)
    !     rei(icol,ilay) = merge(rei_val, 0._wp, iwp(icol,ilay) > 0._wp)

    !   end do
    ! end do
    ! ----------------------------------- icon method ------------------------------------ !
    effective_radius = &
      1.0e6_wp * droplet_scale * (3.0e-9_wp / (4.0_wp * pi * rhoh2o))**(1.0_wp/3.0_wp) 

    reimin = cloud_optics%get_min_radius_ice()
    reimax = cloud_optics%get_max_radius_ice()
    
    relmin = cloud_optics%get_min_radius_liq()
    relmax = cloud_optics%get_max_radius_liq()

    IF (relmax <= relmin .OR. reimax <= reimin) THEN
      CALL stop_on_err('compute_clouds: Droplet minimun size required is bigger than maximum')
    END IF

    !$ACC parallel loop default(none) gang vector collapse(2) async(1)
    do ilay=1,nlay
      do icol=1,ncol
        !
        ! --- Cloud liquid and ice mass: [kg/m2 in cell] --> [g/m2 in cloud]
        !
        ! cld_frc_loc = MAX(EPSILON(1.0_wp),cld_frc(icol,ilay))
        ! ziwp(icol,ilay) = xm_ice(icol,ilay)*1000.0_wp/cld_frc_loc
        ! zlwp(icol,ilay) = xm_liq(icol,ilay)*1000.0_wp/cld_frc_loc

        ! ! Mask which tells cloud optics that this cell is clear
        ! lcldlyr = cld_frc(icol,ilay) > cld_frc_thresh !!!
        ! IF (.NOT. lcldlyr) THEN
        !   ziwp(icol,ilay) = 0.0_wp
        !   zlwp(icol,ilay) = 0.0_wp
        ! END IF

        lwp(icol,ilay) = (p_lev(icol,ilay+1) - p_lev(icol,ilay)) * qc(icol,ilay) / grav * 1e3_wp
        iwp(icol,ilay) = (p_lev(icol,ilay+1) - p_lev(icol,ilay)) * qi(icol,ilay) / grav * 1e3_wp
        !
        ! --- cloud water and ice concentrations [g/m3]
        !
        ziwc = iwp(icol,ilay)/dh_lev(icol,ilay) !!!
        zlwc = lwp(icol,ilay)/dh_lev(icol,ilay)
        !
        ! IF (lcldlyr .AND. (lwp(icol,ilay)+iwp(icol,ilay))>ccwmin) THEN

        IF (lwp(icol,ilay)+iwp(icol,ilay) > ccwmin) THEN

          ! cdnc calculation from icon:
          zprat=(MIN(8._wp,80000._wp/p_lay(icol,ilay)))**2
          zn1 = 0.71_wp*cn1sea + 0.29_wp*cn1lnd
          zn2 = 0.71_wp*cn2sea + 0.29_wp*cn2lnd
          IF (p_lay(icol,ilay).LT.80000._wp) THEN
            zcdnc = 1.e6_wp*(zn1+(zn2-zn1)*(EXP(1._wp-zprat)))
          ELSE
            zcdnc = zn2*1.e6_wp
          END IF
          ! end cdnc calculation from icon:

          zkap = 0.71_wp*zkap_mrtm + 0.29_wp*zkap_cont
          ! todo: extract land mask somehow
          ! IF ( laland(icol) .AND. .NOT.laglac(icol) ) zkap = zkap_cont
          
          re_cryst_scal = MAX(reimin, MIN(reimax,83.8_wp*ziwc**0.216_wp))
          re_drop_scal  = MAX(relmin, MIN(relmax, &
            effective_radius * zkap * (zlwc / zcdnc)**(1.0_wp/3.0_wp) ))

          rei (icol,ilay) = re_cryst_scal
          rel (icol,ilay) = re_drop_scal
        ELSE
          rei (icol,ilay) = reimin
          rel (icol,ilay) = relmin
        END IF
      END DO
    END DO
    ! print *, sum(zcdnc,1)/size(zcdnc,1)

    !$acc exit data delete(cloud_mask)
    !$omp target exit data map(release:cloud_mask)
   
  end subroutine compute_clouds
  !
  ! --------------------------------------------------------------------------------------
  !
  subroutine compute_aerosols
    real(wp), dimension(ncol,nlay) :: vmr_h2o ! h2o vmr
    logical :: is_sulfate, is_dust, is_even_column 
    ! 
    ! Variable and memory allocation 
    !
    if(is_sw) then
      allocate(ty_optical_props_2str::aerosols)
    else
      allocate(ty_optical_props_1scl::aerosols)
    end if
    call stop_on_err(aerosols%init(k_dist%get_band_lims_wavenumber()))
    select type(aerosols)
      class is (ty_optical_props_1scl)
        call stop_on_err(aerosols%alloc_1scl(ncol, nlay))
        !$acc        enter data copyin(aerosols) create(aerosols%tau)
        !$omp target enter data map              (alloc:aerosols%tau)
      class is (ty_optical_props_2str)
        call stop_on_err(aerosols%alloc_2str(ncol, nlay))
        !$acc        enter data copyin(aerosols) create(aerosols%tau, aerosols%ssa, aerosols%g)
        !$omp target enter data map              (alloc:aerosols%tau, aerosols%ssa, aerosols%g)
      class default
        call stop_on_err("my_rte_rrtmgp_allsky: Don't recognize the kind of optical properties ")
    end select
    !
    ! Derive relative humidity from profile
    !   Keep vmr_h2o on the GPU
    !
    !$acc        data create(   vmr_h2o)
    !$omp target data map(alloc:vmr_h2o)
    call stop_on_err(gas_concs%get_vmr("h2o",vmr_h2o))
    !
    ! Aerosol properties 
    ! 
    allocate(aero_type(ncol,nlay), aero_size(ncol,nlay), &
             aero_mass(ncol,nlay), relhum   (ncol,nlay))
    !$acc        enter data create(   aero_type, aero_size, aero_mass, relhum)
    !$omp target enter data map(alloc:aero_type, aero_size, aero_mass, relhum)
    call get_relhum(ncol, nlay, p_lay, t_lay, vmr_h2o, relhum)
    if (do_aerosols) then
      print *, "Clipping relhum to 0, 1 for aerosol computations"
      relhum = min(max(relhum, 0._wp), 1._wp)
    end if
    !$acc end data
    !$omp end target data

    ! Restrict sulfate aerosols to lower stratosphere (> 50 hPa = 50*100 Pa; < 100 hPa = 100*100 Pa)
    !   and dust aerosols to the lower troposphere (> 700 hPa; < 900 hPa), and
    !   put them in 1/2 of the columns
    !
    !
    !$acc                         parallel loop    collapse(2) copyin(p_lay) 
    !$omp target teams distribute parallel do simd collapse(2) map(to:p_lay) 
    do ilay=1,nlay
      do icol=1,ncol
        is_sulfate = (p_lay(icol,ilay) >  50._wp * 100._wp .and. & 
                      p_lay(icol,ilay) < 100._wp * 100._wp)
        is_dust    = (p_lay(icol,ilay) > 700._wp * 100._wp .and. & 
                      p_lay(icol,ilay) < 900._wp * 100._wp)
        is_even_column = mod(icol, 2) /= 0
        if      (is_even_column .and. is_sulfate) then 
          aero_type(icol,ilay) = merra_aero_sulf
          aero_size(icol,ilay) = 0.2_wp
          aero_mass(icol,ilay) = 1.e-6_wp
        else if(is_even_column .and. is_dust) then 
          ! Dust aerosol
          aero_type(icol,ilay) = merra_aero_dust
          aero_size(icol,ilay) = 0.5_wp
          aero_mass(icol,ilay) = 3.e-5_wp
        else
          aero_type(icol,ilay) = 0
          aero_size(icol,ilay) = 0._wp
          aero_mass(icol,ilay) = 0._wp
        end if
      end do
    end do

    end subroutine compute_aerosols
  ! --------------------------------------------------------------------------------------
  !
  ! Calculate layer relative humidity for aerosol optics calculations
  !
  subroutine get_relhum(ncol, nlay, p_lay, t_lay, vmr_h2o, relhum)
    use mo_rte_kind,           only: wp
    use mo_gas_optics_constants,   only: m_h2o, m_dry

    integer,  intent(in) :: ncol, nlay
    real(wp), intent(in) :: p_lay(ncol,nlay)    ! layer pressure (Pa)
    real(wp), intent(in) :: t_lay(ncol,nlay)    ! layer temperature (K)
    real(wp), intent(in) :: vmr_h2o(ncol,nlay)  ! water volume mixing ratio

    real(wp), intent(inout) :: relhum(ncol,nlay) ! relative humidity (fraction, 0-1)

    ! Local variables 
    integer :: i, k

    real(wp) :: mmr_h2o          ! water mass mixing ratio
    real(wp) :: q_lay            ! water specific humidity
    real(wp) :: q_lay_min, q_tmp, es_tmp
    real(wp) :: mwd, t_ref, rh

    ! Set constants
    mwd       = m_h2o/m_dry      ! ratio of water to dry air molecular weights
    t_ref     = 273.16_wp        ! reference temperature (K)
    q_lay_min = 1.e-7_wp         ! minimum water mass mixing ratio
    ! -------------------

    ! Derive layer virtual temperature
    !$acc                         parallel loop    collapse(2) copyin(p_lay, vmr_h2o, t_lay) copyout( relhum)
    !$omp target teams distribute parallel do simd collapse(2) map(to:p_lay, vmr_h2o, t_lay) map(from:relhum) 
    do i = 1, ncol 
       do k = 1, nlay
          ! Convert h2o vmr to mmr
          mmr_h2o = vmr_h2o(i,k) * mwd
          q_lay = mmr_h2o / (1 + mmr_h2o)
          q_tmp = max(q_lay_min, q_lay)
          es_tmp = exp( (17.67_wp * (t_lay(i,k)-t_ref)) / (t_lay(i,k)-29.65_wp) )
          rh = (0.263_wp * p_lay(i,k) * q_tmp) / es_tmp
          ! Convert rh from percent to fraction
          relhum(i,k) = 0.01_wp * rh
       enddo
    enddo
  end subroutine get_relhum
  !--------------------------------------------------------------------------------------------------------------------
  subroutine write_fluxes 
    use netcdf
    use mo_simple_netcdf, only: write_field
    integer :: ncid, i, col_dim, lay_dim, lev_dim, varid, ngpt_dim
    real(wp) :: vmr(ncol, nlay)
    character(len=3) :: flux_prefix
    !
    ! Write fluxes - make this optional? 
    !

    !
    ! Define dimensions 
    !
    if(nf90_create(trim(output_file),  NF90_CLOBBER, ncid) /= NF90_NOERR) &
      call stop_on_err("rrtmgp_allsky: can't create file " // trim(output_file))

    if(nf90_def_dim(ncid, "col", ncol, col_dim) /= NF90_NOERR) &
      call stop_on_err("rrtmgp_allsky: can't define col dimension")
    if(nf90_def_dim(ncid, "lay", nlay, lay_dim) /= NF90_NOERR) &
      call stop_on_err("rrtmgp_allsky: can't define lay dimension")
    if(nf90_def_dim(ncid, "lev", nlay+1, lev_dim) /= NF90_NOERR) &
      call stop_on_err("rrtmgp_allsky: can't define lev dimension")
    
    !
    ! Define variables 
    !
    ! State 
    !
    call create_var("p_lev", ncid, [col_dim, lev_dim])
    call create_var("t_lev", ncid, [col_dim, lev_dim])
    call create_var("p_lay", ncid, [col_dim, lay_dim])
    call create_var("t_lay", ncid, [col_dim, lay_dim])
    call create_var("h2o",   ncid, [col_dim, lay_dim])
    call create_var("o3",    ncid, [col_dim, lay_dim])

    ! All the gases except h2o, o3 - write as attributes? Or not bother? 

    if(do_clouds) then 
      call create_var("lwp", ncid, [col_dim, lay_dim])
      call create_var("iwp", ncid, [col_dim, lay_dim])
      call create_var("rel", ncid, [col_dim, lay_dim])
      call create_var("rei", ncid, [col_dim, lay_dim])
    end if 
    if(do_aerosols) then 
      if(nf90_def_var(ncid, "aero_type", NF90_SHORT, [col_dim, lay_dim], varid) /= NF90_NOERR) &
        call stop_on_err("create_var: can't define variable aero_type")
      call create_var("aero_size", ncid, [col_dim, lay_dim])
      call create_var("aero_mass", ncid, [col_dim, lay_dim])
   end if 
    !
    ! Fluxes - definitions 
    !
    if(is_sw) then 
      flux_prefix = "sw_"
      call create_var(flux_prefix // "flux_dir", ncid, [col_dim, lev_dim])
      ! print *, ngpt
      ! if(nf90_def_dim(ncid, "gpt", ngpt, ngpt_dim) /= NF90_NOERR) &
      !   call stop_on_err("rrtmgp_allsky: can't define ngpt dimension")
      ! call create_var("toa_flux",    ncid, [col_dim, ngpt_dim])
    else
      flux_prefix = "lw_"  
    end if 
    call create_var(flux_prefix // "flux_up", ncid, [col_dim, lev_dim])
    call create_var(flux_prefix //"flux_dn", ncid, [col_dim, lev_dim])
    if(nf90_enddef(ncid) /= NF90_NOERR) &
      call stop_on_err("rrtmgp_allsky: can't end file definition??")

    !
    ! Write variables
    !
    ! State - writing 
    !$acc        update host(p_lev, t_lev, p_lay, t_lay)
    !$omp target update from(p_lev, t_lev, p_lay, t_lay)
    call stop_on_err(write_field(ncid, "p_lev",  p_lev))
    call stop_on_err(write_field(ncid, "t_lev",  t_lev))
    call stop_on_err(write_field(ncid, "p_lay",  p_lay))
    call stop_on_err(write_field(ncid, "t_lay",  t_lay))
    ! Array vmr is on the host, not the device, but is copied-out
    call stop_on_err(gas_concs%get_vmr("h2o", vmr))
    call stop_on_err(write_field(ncid, "h2o",    vmr))
    call stop_on_err(gas_concs%get_vmr("o3",  vmr))
    call stop_on_err(write_field(ncid, "o3",     vmr))

    if(do_clouds) then 
      !$acc        update host(lwp, iwp, rel, rei)
      !$omp target update from(lwp, iwp, rel, rei)
      call stop_on_err(write_field(ncid, "lwp",  lwp))
      call stop_on_err(write_field(ncid, "iwp",  iwp))
      call stop_on_err(write_field(ncid, "rel",  rel))
      call stop_on_err(write_field(ncid, "rei",  rei))
    end if 

    if(do_aerosols) then 
      !$acc        update host(aero_size, aero_mass, aero_type)
      !$omp target update from(aero_size, aero_mass, aero_type)
      call stop_on_err(write_field(ncid, "aero_size",  aero_size))
      call stop_on_err(write_field(ncid, "aero_mass",  aero_mass))
      call stop_on_err(write_field(ncid, "aero_type",  aero_type))
    end if 

    ! if(is_sw) then
    !   call stop_on_err(write_field(ncid, "toa_flux",  toa_flux))
    ! end if
    ! Fluxes - writing 
    !$acc        update host(flux_up, flux_dn)
    !$omp target update from(flux_up, flux_dn)
    call stop_on_err(write_field(ncid, flux_prefix // "flux_up",  flux_up))
    call stop_on_err(write_field(ncid, flux_prefix // "flux_dn",  flux_dn))
    if(.not. is_lw) then 
      !$acc        update host(flux_dir)
      !$omp target update from(flux_dir)
      call stop_on_err(write_field(ncid, flux_prefix // "flux_dir",  flux_dir))
    end if 

    ! Close netCDF 
    if(nf90_close(ncid) /= NF90_NOERR) call stop_on_err("rrtmgp_allsky: error closing file??")
  end subroutine write_fluxes
  ! ---------------------------------------------------------
  subroutine create_var(name, ncid, dim_ids)
    use netcdf
    character(len=*),      intent(in) :: name 
    integer,               intent(in) :: ncid
    integer, dimension(:), intent(in) :: dim_ids 

    integer :: varid

    if(nf90_def_var(ncid, trim(name), NF90_DOUBLE, dim_ids, varid) /= NF90_NOERR) &
      call stop_on_err("create_var: can't define " // trim(name) // " variable")
  end subroutine create_var
  ! ---------------------------------------------------------
  subroutine get_gas_vmr(ncol, nlay, mmr, gas_name, vmr)
  !from mmr to vmr (see https://github.com/E3SM-Project/ACME-ECP/blob/master/components/cam/src/physics/rrtmgp/radiation.F90)

      integer, intent(in ) :: ncol, nlay
      character(len=*), intent(in) :: gas_name
      real(wp), intent(in) :: mmr(ncol,nlay)
      real(wp), intent(out) :: vmr(ncol,nlay)

      ! vmr = 1.0_wp

      ! Gases and molecular weights. Note that we do NOT have CFCs yet (I think
      ! this is coming soon in RRTMGP). RRTMGP also allows for absorption due to
      ! CO and N2, which RRTMG did not have.
      character(len=3), dimension(8) :: gas_species = (/ &
        'H2O', 'CO2', 'O3 ', 'N2O', &
        'CO ', 'CH4', 'O2 ', 'N2 ' &
      /)
      real(wp), dimension(8) :: mol_weight_gas = (/ &
        18.01528, 44.0095, 47.9982, 44.0128, &
        28.0101, 16.04246, 31.998, 28.0134 &
      /)  ! g/mol

      ! Molar weight of air
      real(wp), parameter :: mol_weight_air = 28.97  ! g/mol
                                      
      ! Defaults for gases that are not available (TODO: is this still accurate?)
      real(wp), parameter :: co_vmr = 1.0e-7_wp
      real(wp), parameter :: n2_vmr = 0.7906_wp

      ! Loop indices
      integer :: igas

      ! Name of routine
      character(len=32) :: subname = 'get_gas_vmr'

      ! Get index into gas names we define above that we know the molecular 
      ! weights for; if this gas is not in list of gases we know about, skip
      igas = string_loc_in_array(gas_name, gas_species)
      !igas = 1
      if (igas <= 0) then
        !call endrun('Gas name ' // trim(gas_name) // ' not recognized.')
        call stop_on_err('Gas name ' // trim(gas_name) // ' not recognized.')
      end if
        
      ! initialize
      vmr(:,:) = 0._wp

      select case(trim(gas_species(igas)))

        case('CO')

           ! CO not available, use default
           vmr(1:ncol,1:nlay) = co_vmr

        case('N2')

           ! N2 not available, use default
           vmr(1:ncol,1:nlay) = n2_vmr

        case('H2O')
           ! Convert to volume mixing ratio by multiplying by the ratio of
           ! molecular weight of dry air to molecular weight of gas. Note that
           ! first specific humidity (held in the mass_mix_ratio array read
           ! from rad_constituents) is converted to an actual mass mixing
           ! ratio.
           vmr(1:ncol,1:nlay) = mmr(1:ncol,1:nlay) / ( &
              1._wp - mmr(1:ncol,1:nlay) &
           )  * mol_weight_air / mol_weight_gas(igas)

        case DEFAULT

           ! Convert to volume mixing ratio by multiplying by the ratio of
           ! molecular weight of dry air to molecular weight of gas
           vmr(1:ncol,1:nlay) = mmr(1:ncol,1:nlay) &
                              * mol_weight_air / mol_weight_gas(igas)

      end select

   end subroutine get_gas_vmr
end program my_rte_rrtmgp_allsky
