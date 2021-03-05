#!/usr/bin/bash

DICE_EXE=/projects/jasm3285/apps/pull_requests/Dice/Dice

# Generate all Dice input Files
python gen_input.py

# Run all Dice jobs
for i in {0..9}
do
    mpirun -np ${OMP_NUM_THREADS} $DICE_EXE input_${i}.dat > output_${i}.dat
done