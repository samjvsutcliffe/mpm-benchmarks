#! /bin/bash
python3 pullout_gen.py
echo "Removing files"
rm pullout/results/pullout/*
echo "Running code"
mpm -p 8 -i ./pullout/input_file.json -f ./ | tee out.txt
