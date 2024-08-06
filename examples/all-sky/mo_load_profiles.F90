module mo_load_profiles
  use mo_rte_kind,      only: wp
  ! use mo_optical_props, only: ty_optical_props,      &
  !                             ty_optical_props_arry, &
  !                             ty_optical_props_1scl, &
  !                             ty_optical_props_2str, &
  !                             ty_optical_props_nstr
  ! use mo_cloud_optics_rrtmgp, & 
  !                       only: ty_cloud_optics_rrtmgp
  use mo_simple_netcdf, only: read_field, read_string, var_exists, get_dim_size, &
                              write_field, create_dim, create_var
  use netcdf

  implicit none
  private
  public :: load_profiles, load_vert_grid_info
  ! ----------------------------------------------------------------------------------

contains
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! read profiles from netcdf
  !
  subroutine load_profiles(cs_mli_file, ncol, nlay, hyai, hybi, hyam, hybm, p0, p_lay, t_lay, p_lev, t_lev, q_lay, qc_lay, qi_lay, o3, n2o, ch4, lwup)
    ! -----------------
    character(len=256),                intent(in)  :: cs_mli_file
    integer,                           intent(in)  :: ncol, nlay
    real(wp), dimension(nlay+1),       intent(in)  :: hyai, hybi
    real(wp), dimension(nlay),         intent(in)  :: hyam, hybm
    real(wp),                          intent(in)  :: p0
    real(wp), dimension(ncol, nlay  ), intent(out) :: p_lay, t_lay, q_lay, qc_lay, qi_lay, o3, n2o, ch4
    real(wp), dimension(ncol, nlay+1), intent(out) :: p_lev, t_lev
    real(wp), dimension(ncol),         intent(out) :: lwup

    ! Local variables
    integer :: ncid, i
    real(wp), dimension(ncol) :: ps

    ! real(wp), dimension(:,:), allocatable                :: q_lay
    ! -----------------
    ! Open cloud optical property coefficient file
    if(nf90_open(trim(cs_mli_file), NF90_NOWRITE, ncid) /= NF90_NOERR) &
       call stop_on_err("load_cld_lutcoeff(): can't open file " // trim(cs_mli_file))
    ! if(nf90_open(trim(cld_coeff_file), NF90_NOWRITE, ncid) /= NF90_NOERR) &
    !    call stop_on_err("load_cld_lutcoeff(): can't open file " // trim(cld_coeff_file))

    ps = read_field(ncid, 'state_ps', 384)
    ! p_lay = hyam * p0 + hybm * ps
    ! p_lev = hyai * p0 + hybi * ps
    call compute_outer_product(ps, hybm, ncol, nlay, p_lay)
    call compute_outer_product(ps, hybi, ncol, nlay+1, p_lev)
    do i=1, ncol
      p_lay(i,:) = p_lay(i,:) + hyam * p0
      p_lev(i,:) = p_lev(i,:) + hyai * p0
    end do

    t_lay = read_field(ncid, 'state_t', 384, 60)
    q_lay = read_field(ncid, 'state_q0001', 384, 60)
    qc_lay = read_field(ncid, 'state_q0002', 384, 60)
    qi_lay = read_field(ncid, 'state_q0003', 384, 60)
    o3 = read_field(ncid, 'pbuf_ozone', 384, 60)
    n2o = read_field(ncid, 'pbuf_N2O', 384, 60)
    ch4 = read_field(ncid, 'pbuf_CH4', 384, 60)
    lwup = read_field(ncid, 'cam_in_LWUP', 384)

    ncid = nf90_close(ncid)
    ! call stop_on_err(cloud_spec%load(band_lims_wvn,                      &
    !                                  radliq_lwr, radliq_upr,             &
    !                                  radice_lwr, radice_upr,             &
    !                                  lut_extliq, lut_ssaliq, lut_asyliq, &
    !                                  lut_extice, lut_ssaice, lut_asyice))
  end subroutine load_profiles
  ! ---------------------------------
  subroutine load_vert_grid_info(cs_grid_info_file, ncol, nlay, hyai, hybi, hyam, hybm, p0)
    ! -----------------
    character(len=256),                intent(in)  :: cs_grid_info_file
    integer,                           intent(in ) :: ncol, nlay
    real(wp), dimension(nlay+1), intent(out)         :: hyai, hybi
    real(wp), dimension(nlay), intent(out)       :: hyam, hybm
    real(wp), intent(out)                          :: p0

    ! Local variables
    integer :: ncid

    if(nf90_open(trim(cs_grid_info_file), NF90_NOWRITE, ncid) /= NF90_NOERR) &
       call stop_on_err("load_vert_grid_info(): can't open file " // trim(cs_grid_info_file))

    hyai = read_field(ncid, 'hyai', nlay+1)
    hybi = read_field(ncid, 'hybi', nlay+1)
    hyam = read_field(ncid, 'hyam', nlay)
    hybm = read_field(ncid, 'hybm', nlay)
    p0 = read_field(ncid, 'P0')


    ncid = nf90_close(ncid)
    ! call stop_on_err(cloud_spec%load(band_lims_wvn,                      &
    !                                  radliq_lwr, radliq_upr,             &
    !                                  radice_lwr, radice_upr,             &
    !                                  lut_extliq, lut_ssaliq, lut_asyliq, &
    !                                  lut_extice, lut_ssaice, lut_asyice))
  end subroutine load_vert_grid_info
  ! -----------------------------------------------------------------------------------
    subroutine stop_on_err(msg)
      !
      ! Print error message and stop
      !
      use iso_fortran_env, only : error_unit
      character(len=*), intent(in) :: msg
      if(len_trim(msg) > 0) then
        write (error_unit,*) trim(msg)
        error stop 1
      end if
    end subroutine

    subroutine compute_outer_product(a, b, n, m, C)
        integer, intent(in) :: n, m
        real(wp), intent(in) :: a(n), b(m)
        real(wp), intent(out) :: C(n, m)
        integer :: i, j

        ! Compute the outer product
        do i = 1, n
            do j = 1, m
              ! print *,C(i,j)
              ! print *,a(i)
              ! print *,b(j)
              ! print *,i,j
              C(i,j) = a(i) * b(j)
            end do
        end do
    end subroutine compute_outer_product

end module
