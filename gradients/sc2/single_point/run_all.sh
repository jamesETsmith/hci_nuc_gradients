#!/bin/bash

# Setup directories
mkdir -p _figures _data _logs

# Generate data
python sc2_run_all.py
python plot_grad.py

# Check for converged in MCSCF
echo "Calculations that converged successfully"
echo "========================================"
grep "CASSCF converged" _logs/*.out

echo ""
echo ""

echo "Calculations that failed to converged"
echo "====================================="
grep "not converged" _logs/*.out

# Clean up
rm -rf *.bkp FCIDUMP *.txt *.dat shci.e
