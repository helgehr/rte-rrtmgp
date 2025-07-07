#!/bin/bash

#SBATCH --partition=develbooster
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --gres=gpu:4
#SBATCH --account=icon-a-ml
#SBATCH --time=02:00:00
#SBATCH --output=slurm_out-comp_tends-%j.out

#if (( $# != 1 )); then
    #echo "Expecting 1 args: <py_script>"
    #exit 1
#fi

echo "Date: $(date)"
echo "Starting py script"

cat ./commands.txt | xargs -I CMD --max-procs=96 bash -c CMD

echo "Exit code of python script $?"
echo "Date: $(date)"
