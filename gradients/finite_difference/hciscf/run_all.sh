#!/usr/bin/bash

# Create data directories
mkdir -p _logs _data _figures

# Generate data for various finite-difference step sizes
for H in 0.1 0.01 0.001 0.0001
do
    echo "Starting h=${H}"
    python n2_casscf_grad.py $H > _logs/h=$H.out
    mv output.dat _data/$H.dat
    echo "Done h=${H}"
done

# Check for convergence problems
grep "not converged" _logs/*.out

# Clean up
rm -rf *.bkp FCIDUMP *.txt *.dat shci.e