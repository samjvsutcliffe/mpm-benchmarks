#! /bin/bash
echo "Removing files"
rm consol/results/consol/*
echo "Running code"
python3 converg_gen.py | tee out.txt
#mpm -p 8 -i ./consol/input_file.json -f ./ | tee out.txt
