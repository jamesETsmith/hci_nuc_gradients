#!/usr/bin/bash

# Set up directories
mkdir -p _molden _logs _data _chk

# Geometry optimization
python ../../opt.py --eps1=5e-05 --ncas=40 --spin=0 --geometry=A 2>_logs/geometric_5e-05.out
