#! /bin/sh
python3 consol_gen.py
echo "Removing files"
rm consol/results/consol/*
export OMP_SCHEDULE="static"
echo "Running code"
mpm -p 8 -i ./consol/input_file.json -f ./
