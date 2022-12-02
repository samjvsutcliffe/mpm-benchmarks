#! /bin/bash
python3 slump_gen.py
echo "Removing files"
rm slump/results/slump/*
echo "Running code"
/mnt/e/Github/mpm/build/mpm -p 8 -i ./slump/input_file.json -f ./
