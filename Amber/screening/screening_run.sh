#!/usr/bin/env bash
#SBATCH -J amber_array
#SBATCH -A your_account
#SBATCH -p your_partition
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=8
#SBATCH -t 24:00:00
#SBATCH --array=1-100%10    # run 100 systems, 10 in parallel
#SBATCH -o logs/%x_%a.out
#SBATCH -e logs/%x_%a.err

module load amber
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

ROOT=$PWD
SYS=$(sed -n "${SLURM_ARRAY_TASK_ID}p" systems.txt)
WORKDIR="${ROOT}/${SYS}"

mkdir -p "${WORKDIR}/build" "${WORKDIR}/mdin" "${WORKDIR}/logs" "${WORKDIR}/run"

# --- write leap.in for this system (could also be pre-copied) ---
cat > "${WORKDIR}/build/leap.in" <<'EOF'
source leaprc.protein.ff14SB
source leaprc.water.tip3p
if ( fileexists("./mcpb.lib") ) { loadOff ./mcpb.lib }
if ( fileexists("./mcpb.frcmod") ) { loadAmberParams ./mcpb.frcmod }
prot = loadPdb ../input.pdb
solvateOct prot TIP3PBOX 10.0
addIonsRand prot Na+ 0
addIonsRand prot Cl- 0
saveAmberParm prot ../run/prmtop ../run/inpcrd
quit
EOF

# --- copy mdin templates (assumed to be in ROOT/mdin_templates) ---
cp -r "${ROOT}/mdin_templates/"* "${WORKDIR}/mdin/" || true

# --- build ---
cd "${WORKDIR}/build"
tleap -f leap.in | tee ../logs/tleap.log

# --- run (pmemd.MPI) ---
cd "${WORKDIR}/run"
srun pmemd.MPI -O -i ../mdin/min_1.in -o min_1.out -p prmtop -c inpcrd    -r min_1.rst7 -ref inpcrd
srun pmemd.MPI -O -i ../mdin/min_2.in -o min_2.out -p prmtop -c min_1.rst7 -r min_2.rst7
srun pmemd.MPI -O -i ../mdin/heat.in   -o heat.out   -p prmtop -c min_2.rst7 -r heat.rst7 -x heat.nc -inf heat.info -ref min_2.rst7
srun pmemd.MPI -O -i ../mdin/eq.in     -o eq.out     -p prmtop -c heat.rst7  -r eq.rst7   -x eq.nc   -inf eq.info   -ref heat.rst7
srun pmemd.MPI -O -i ../mdin/prod.in   -o prod.out   -p prmtop -c eq.rst7    -r prod.rst7 -x prod.nc -inf prod.info
