# 1) Prepare MCPB input (edit selection masks to your coordinating atoms)
mcpb.py -i mcpb.in -s 1    # or the multi-step workflow from the tutorial

# 2) After MCPB finishes, it writes something like: mcpb.lib  mcpb.frcmod  complex_mcpb.pdb
# 3) tLEaP build (bonded model)
cat > build/leap_mcpb.in <<'EOF'
source leaprc.protein.ff14SB
source leaprc.water.tip3p
loadOff  ./mcpb.lib
loadAmberParams ./mcpb.frcmod
mol = loadPdb ../protein_pH7_POCKET001_GD3.pdb
solvateOct mol TIP3PBOX 10.0
addIonsRand mol Cl- 0
addIonsRand mol Na+ 0
saveAmberParm mol ../run/prmtop ../run/inpcrd
quit
EOF

tleap -f build/leap_mcpb.in | tee logs/tleap_mcpb.log
