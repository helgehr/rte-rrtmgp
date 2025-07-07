- activate dev conda environment to, most importantly, set environment variables
- export RRTMGP_ROOT=/path/to/rte-rrtmgp
- export RRTMGP_DATA=/path/to/rrtmgp-data
- export FCINCLUDE=-I${CONDA_PREFIX}/include
# Eventually add "-lstdc++" to LIBS in Makefiles var as: LIBS += -lnetcdff -lnetcdf -lstdc++
