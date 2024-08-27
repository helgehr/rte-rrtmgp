set -eu #x

if (( $# != 4 )); then
	echo "Please provide a climsim file and out_file"
	exit 1
fi

ncol=$1
nlay=$2
mli_file=$3
out_path=$4
py_script="${RRTMGP_ROOT}/examples/all-sky/rad_flux2tend.py"

#ncol="384" #"21600"#
#nlay="60"

if (( $ncol != "384" )) && [[ "$mli_file" =~ "low-res" ]]; then
  echo "Infile seems to be low_res but using $ncol columns now"
  exit 1
fi

if (( $ncol != "21600" )) && [[ "$mli_file" =~ "high-res" ]]; then
  echo "Infile seems to be high_res but using $ncol columns now"
  exit 1
fi

tmp_dir=$(mktemp -d)
#tmp_dir="./testout/"
#echo $tmp_dir
sw_path="${tmp_dir}/my-rrtmgp-allsky-sw_noaero.nc"
lw_path="${tmp_dir}/my-rrtmgp-allsky-lw_noaero.nc"

./my_rrtmgp_allsky $ncol $nlay 1 $sw_path \
		${RRTMGP_DATA}/rrtmgp-gas-sw-g224.nc \
		/p/scratch/icon-a-ml/heuer1/LEAP/ClimSim_low-res/ClimSim_low-res_grid-info.nc \
		$mli_file \
		${RRTMGP_DATA}/rrtmgp-clouds-sw.nc

./my_rrtmgp_allsky $ncol $nlay 1 $lw_path \
		${RRTMGP_DATA}/rrtmgp-gas-lw-g256.nc \
		/p/scratch/icon-a-ml/heuer1/LEAP/ClimSim_low-res/ClimSim_low-res_grid-info.nc \
		$mli_file \
		${RRTMGP_DATA}/rrtmgp-clouds-lw.nc

# # first steps with cdo, but I have to used python anyway for the consecutive differences -> transition to one python script
# cdo expr,'net_flux=sw_flux_up-sw_flux_dn' my-rrtmgp-allsky-sw_noaero.nc my-rrtmgp-allsky-sw_noaero_netflux.nc
# cdo expr,'net_flux=lw_flux_up-lw_flux_dn' my-rrtmgp-allsky-lw_noaero.nc my-rrtmgp-allsky-lw_noaero_netflux.nc
# python xr_var_diff.py my-rrtmgp-allsky-sw_noaero_netflux.nc y my-rrtmgp-allsky-sw_noaero_fluxdiv.nc
# python xr_var_diff.py my-rrtmgp-allsky-lw_noaero_netflux.nc y my-rrtmgp-allsky-lw_noaero_fluxdiv.nc
# cdo add my-rrtmgp-allsky-sw_noaero_fluxdiv.nc my-rrtmgp-allsky-lw_noaero_fluxdiv.nc my-rrtmgp-allsky-noaero_fluxdiv.nc

mamba run -n heuer1_rte_rrtmgp_dev --no-capture-output python -u $py_script $sw_path $lw_path $out_path

rm -rf $tmp_dir
