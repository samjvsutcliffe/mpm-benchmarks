#!/bin/bash

# Request resources:
#SBATCH -c 16     # 1 entire node
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
pip install -e ~/cb-geo/pycbg

source ~/cb-geo/mpm/build/setup-vars.sh
export PATH=~/cb-geo/mpm/build/:$PATH
echo "Removing files"
rm -r ./consol_conv*
echo "Generating files"
python3 converg_gen.py
