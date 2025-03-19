- activate dev conda environment to, most importantly, set environment variables
- load Stages/2024
- export RRTMGP_ROOT=/path/to/rte-rrtmgp
- export RRTMGP_DATA=/path/to/rrtmgp-data
- export FCINCLUDE=-I${CONDA_PREFIX}/include
# Eventually (first worked without, now not anymore), export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/p/project/icon-a-ml/mambaforge/envs/heuer1_rte_rrtmgp_dev/lib"
# Eventually add "-lstdc++" to LIBS in Makefiles var as: LIBS += -lnetcdff -lnetcdf -lstdc++
