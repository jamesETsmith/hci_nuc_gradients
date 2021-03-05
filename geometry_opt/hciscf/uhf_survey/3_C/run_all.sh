#!/bin/bash

# Set up data dirs
mkdir -p _slurm _data _molden _logs _chk

MULT=$1
GEOM=$2

# Run UHF for all three opt strategies
python ../uhf.py ${MULT} ${GEOM} diis
python ../uhf.py ${MULT} ${GEOM} adiis
python ../uhf.py ${MULT} ${GEOM} newton

