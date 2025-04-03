
gm pdb2gmx -f protein.pdb -o processed.gro -ff 15 -water spce
gmx editconf -f processed.gro -o newbox.gro -c -d 1.0 -bt cubic
gmx solvate -cp newbox.gro -cs spc216.gro -o solved.gro -p topol.top
gmx grompp -f ions.mdp
