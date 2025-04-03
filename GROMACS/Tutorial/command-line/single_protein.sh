
gm pdb2gmx -f protein.pdb -o processed.gro -ff 15 -water spce
# It is advise to take a  box, but a cube is fine for the start.
gmx editconf -f processed.gro -o newbox.gro -c -d 1.0 -bt cubic
gmx solvate -cp newbox.gro -cs spc216.gro -o solved.gro -p topol.top
gmx grompp -f ions.mdp -c solved.gro -p topol.top -ions.tpr

# Add the ions name after -pname (for the positiv ion) and -nname (for the negativ ion). -pq to specify the chare of the positive ion and -nq for the charge of the negativ ions.
gmx genion -s ions.tpr -o solved_ions.gro -p topol.top -pname NA -nname CL -neutral
# Build the file for the Energy minimization the minim.mdp should be written by the author.
gmx grompp -f minim.mdp -c solved_ions.gro -p topol.top -o em.tpr
gmx mdrun -deffnm em
gmx energy -f em.edr -o potential.xvg # To produce the energy minimization graph.

# Calculation of the equilibration respective to the temperature, the nvt.mdp file shoud be written by the author 
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -deffnm nvt
gmx energy -f nvt.edr -o temperature.xvg # To produce the graph.

# Calculation of the equilibration respective to the pressure, the npt.mdp file shoud be written by the author 
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -npt.tpr
gmx mdrun -deffnm npt
gmx energy -f npt.edr -o pressure.xvg
