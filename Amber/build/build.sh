#!/usr/bin/env bash
set -euo pipefail
mkdir -p ../run ../logs
tleap -f leap.in | tee ../logs/tleap.log
echo "Built prmtop/inpcrd in run/"
