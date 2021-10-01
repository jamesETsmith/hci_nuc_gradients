#!/bin/bash

# 
# for D in DIIS ADIIS FN_DIIS FN_ADIIS NEWTON
for D in NEWTON
do
	cd $D
	# cd 1_A && ../../submit.sh $D 1 A && cd ..
	# cd 1_B && ../../submit.sh $D 1 B && cd ..
	cd 1_C && ../../submit.sh $D 1 C && cd ..
	# cd 3_A && ../../submit.sh $D 3 A && cd ..
	cd 3_B && ../../submit.sh $D 3 B && cd ..
	cd 3_C && ../../submit.sh $D 3 C && cd ..
	cd ..
done
