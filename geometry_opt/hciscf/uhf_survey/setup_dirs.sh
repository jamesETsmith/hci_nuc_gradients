#!/bin/bash

for D in DIIS ADIIS FN_DIIS FN_ADIIS NEWTON
do 
	mkdir -p $D
	cd $D
	mkdir -p 1_A 1_B 1_C 3_A 3_B 3_C	
	cd ..
done
