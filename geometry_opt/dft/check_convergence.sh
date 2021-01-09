#!/bin/bash

echo "Converged"
echo "========="
grep --color -m 1 "Frequencies --" ./*/_opt.log

echo ""
echo "Not converged"
echo "============="
grep -L "Frequencies --" ./*/_opt.log