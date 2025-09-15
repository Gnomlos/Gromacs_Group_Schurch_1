#!/usr/bin/env bash
#SBATCH -J amber_gpu
#SBATCH -A your_account
#SBATCH -p gpu
#SBATCH -N 1
#SBATCH --gres=gpu:1
#SBATCH -t 12:00:00
#SBATCH -o logs/amber.out
#SBATCH -e logs/amber.err

module load amber/cuda   # site-specific
cd "$(dirname "$0")"

PMEMD=pmemd.cuda

$PMEMD -O -i ../mdin/min_1.in -o min_1.out -p prmtop -c inpcrd -r min_1.rst7 -ref inpcrd
$PMEMD -O -i ../mdin/min_2.in -o min_2.out -p prmtop -c min_1.rst7 -r min_2.rst7
$PMEMD -O -i ../mdin/heat.in   -o heat.out   -p prmtop -c min_2.rst7 -r heat.rst7 -x heat.nc -inf heat.info -ref min_2.rst7
$PMEMD -O -i ../mdin/eq.in     -o eq.out     -p prmtop -c heat.rst7  -r eq.rst7   -x eq.nc   -inf eq.info   -ref heat.rst7
$PMEMD -O -i ../mdin/prod.in   -o prod.out   -p prmtop -c eq.rst7    -r prod.rst7 -x prod.nc -inf prod.info
