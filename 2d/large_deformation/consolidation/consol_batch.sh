#! /bin/bash
python3 consol_gen.py
echo "Removing files"
rm consol/results/consol/*
echo "Running code"
~/cb-geo/mpm/build/mpm -p 8 -i ./consol/input_file.json -f ./ | tee out.txt
