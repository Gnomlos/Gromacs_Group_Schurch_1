gmx pdb2gmx -f protein.pdb -o processed.gro -ff OPLS-AA -water spce
# It is advise to take a  box, but a cube is fine for the start.
gmx editconf -f processed.gro -o newbox.gro -c -d 1.0 -bt cubic
gmx solvate -cp newbox.gro -cs spc216.gro -o solved.gro -p topol.top
gmx grompp -f ions.mdp -c solved.gro -p topol.top -ions.tpr

# Add the ions name after -pname (for the positiv ion) and -nname (for the negativ ion). -pq to specify the chare of the positive ion and -nq for the charge of the negativ ions.
gmx genion -s ions.tpr -o solved_ions.gro -p topol.top -pname NA -nname CL -neutral
# Build the file for the Energy minimization the minim.mdp should be written by the author.
gmx grompp -f minim.mdp -c solved_ions.gro -p topol.top -o em.tpr
gmx mdrun -deffnm em
echo Pot | gmx energy -f em.edr -o potential.xvg # To produce the energy minimization graph. # To draw the graph

# Calculation of the equilibration respective to the temperature, the nvt.mdp file shoud be written by the author 
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -deffnm nvt
echo Temp | gmx energy -f nvt.edr -o temperature.xvg # To produce the graph.

# Calculation of the equilibration respective to the pressure, the npt.mdp file shoud be written by the author 
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -npt.tpr
gmx mdrun -deffnm npt
echo Pres | gmx energy -f npt.edr -o pressure.xvg

#Calculating density
echo Dens | gmx energy -f npt.edr -o density.xvg

#Simulation of the molecule after equilibration
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0_1.tpr
gmx mdrun -deffnm md_0_1

# Analysis
echo Prot Sys | gmx trjconv -s md_0_1.tpr -f md_0_1.xtc -o md_0_1_noPBC.xtc -pbc mol -center
echo Back Back | gmx rms -s md_0_1.tpr -f md_0_1_noPBC.xtc -o rmsd.xvg -tu ns
echo Prot | gmx gyrate -s md_0_1.tpr -f md_0_1_noPBC.xtc -o gyrate.xvg






