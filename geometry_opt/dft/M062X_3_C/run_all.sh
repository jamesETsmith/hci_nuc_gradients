#!/bin/bash

# User settings
file_root=${PWD##*/}

# Run the single point
g16 sp.gjf _sp.log

# Run the geometry optimization
g16 opt.gjf _opt.log

# Run natural orbital analysis 
g16 natorb.gjf _natorb.log

# Convert to fchk
formchk ${file_root}_NO.chk ${file_root}.fchk

# Collect the natural orbital occupation numbers
python ../parse_for_noons.py _natorb.log

# Generate cube files
mkdir -p _cubes
echo "Writing cube files for ${file_root}"
cubegen 24 MO=118 ${file_root}.fchk _cubes/${file_root}_118.cub
cubegen 24 MO=119 ${file_root}.fchk _cubes/${file_root}_119.cub
cubegen 24 MO=120 ${file_root}.fchk _cubes/${file_root}_120.cub
cubegen 24 MO=121 ${file_root}.fchk _cubes/${file_root}_121.cub

# Clean up
rm -rf Gau-*