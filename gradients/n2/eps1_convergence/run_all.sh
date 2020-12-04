#!/bin/bash

# Setup directories
mkdir -p _figures _data _logs

# Generate data
python vhciscf.py
python vhciscf_aa.py
python hciscf.py
python hciscf_aa.py

# Plot results
python plot_conv.py


# Check for converged in MCSCF
echo ""
echo ""
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
