#!/bin/bash

# Request resources:
#SBATCH -n 4     # 1 entire node
#SBATCH -N 1     # 1 entire node
#SBATCH -c 1     # 1 entire node
#SBATCH --time=06:00:0  # 6 hours (hours:minutes:seconds)
#SBATCH --mem=2G      # 1 GB RAM
#SBATCH -p shared


echo "Running code"

module load gcc
module load mkl
module load boost
module load eigen
module load openmpi
module load python
pip install pycbg
source ~/cb-geo/mpm/build/setup-vars.sh

echo "Generating files"
python3 slump_gen.py
echo "Removing files"
rm ./slump/results/slump/*
echo "Running code"
mpirun mpm -p 1 -f ./ -i ./slump/input_file.json
