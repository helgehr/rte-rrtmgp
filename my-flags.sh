#export FC="gfortran"
#export FCFLAGS="-ffree-line-length-none -m64 -std=f2008 -march=native -fbounds-check -fmodule-private -fimplicit-none -finit-real=nan -g -DRTE_USE_CBOOL"
export RRTMGP_ROOT="/p/project/icon-a-ml/heuer1/RadiationSchemes/rte-rrtmgp"
export RRTMGP_DATA="/p/project1/icon-a-ml/heuer1/RadiationSchemes/rrtmgp-data"
export FCINCLUDE=-I${CONDA_PREFIX}/include
