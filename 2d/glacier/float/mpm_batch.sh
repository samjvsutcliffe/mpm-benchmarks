#!/bin/bash

# Request resources:
#SBATCH -c 32     # 1 entire node
#SBATCH --time=10:00:0  # 6 hours (hours:minutes:seconds)
#SBATCH --mem=3G      # 1 GB RAM
#SBATCH -p shared


echo "Running code"

module load gcc
#module load mkl
module load openblas
module load boost
module load eigen
module load python
pip install pycbg

source ~/cb-geo/mpm/build/setup-vars.sh
echo "Generating files"
python3 slump_gen.py
echo "Removing files"
rm ./notch/results/notch/*
echo "Running code"
~/cb-geo/mpm/build/mpm -f ./ -i ./notch/input_file.json
