#!/bin/bash

# Run UHF and generate NOs (Step #1-#3 in the paper)
python uhf.py

# Use the UNO orbitals to run a large and loose HCISCF
python vhciscf.py

