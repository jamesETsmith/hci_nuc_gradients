#!/usr/bin/bash

# Set up directories
mkdir -p _molden _logs _data _chk

# Geometry optimization
EPS1=2e-4
NCAS=40

python ../opt.py ${EPS1} ${NCAS} 1> _logs/opt_${EPS1}.out 2>_logs/geometric_${EPS1}.out

