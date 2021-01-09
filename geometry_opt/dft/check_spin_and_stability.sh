#!/bin/bash

# Checking SP calculations
for f in $(find _sp_logs/ -name "*.log"); do echo $f && grep --color --  "<S\*\*2>" $f | tail -2; done

for f in $(find _sp_logs/ -name "*.log"); do echo $f && grep -A 2 --color --  "spin densiti" $f | tail -2; done


#