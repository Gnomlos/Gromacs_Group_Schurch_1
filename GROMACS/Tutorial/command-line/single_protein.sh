
gm pdb2gmx -f protein.pdb -o processed.gro -ff 15 -water spce
gmx editconf -f processed.gro -o newbox.gro -c -d 1.0 -bt cubic
gmx solvate -cp newbox.gro -cs spc216.gro -o solved.gro -p topol.top
gmx grompp -f ions.mdp -c solved.gro -p topol.top -ions.tpr
# Add the ions name after -pname (for the positiv ion) and -nname (for the negativ ion). -pq to specify the chare of the positive ion and -nq for the charge of the negativ ions.
gmx genion -s ions.tpr -o solved_ions.gro -p topol.top -pname NA -nname CL -neutral

