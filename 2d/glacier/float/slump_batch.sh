#! /bin/bash
python3 slump_gen.py
echo "Removing files"
rm notch/results/notch/*
echo "Running code"
export OMP_SCHEDULE="static"
#/mnt/d/Github/mpm/build/mpm -i ./slump/input_file.json -f ./
mpm -p 4 -i ./notch/input_file.json -f ./
#| tee out.txt
