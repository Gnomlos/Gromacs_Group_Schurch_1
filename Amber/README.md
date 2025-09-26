# AMBER: A Molecular dynamic simulation software
Amber, much like GROMACS is one of the most used md software. The main difference with GROMACS is however its **flexibility**. Over time multiple settings and tools for different case were devlopped. In our case, their is two main tools
wich are very intersesting: MCBP.py and MMPBSA.py. Both of wich will serve us to simulate the Gd ions and to calculate the bonding energy of this ion to the protein (possibly molecule). 
> For the installation of Amber, please consult the instruction in the *Amber_source_file* folder.
## A quick introduction into molecular dynamic
Before starting using Amber it is preferable to understand what exactly we are doing. This will helps us to better interpret our results and check for any problem as well as correcting problems. **Molecular dynamic simulation** is *simply*
the simulation of a chemical environement using a computer. However contrary to quantum chemistry computin (like **ORCA**), it is based on classical physics equation and not on quantum theory. The advantage of this is that it allows for a 
simpler, lighter software which gives overall good results; it is however a bit less precise. 

Ok, but how do we simulate our experience? Simply buy telling the software what are the conditions of our experiment **in a precise order**. The order of steps is important and it is actually where the major errors can come from; so lets decompose the different steps and explain what each do. 

### System Preparation
This step is used to setup the environement for our molecule of study. So image you had a glass cube and you wanted to study a protein. First you would put the protein in the cube then fill it with the desire solvent (in most of the case water) Then defining your pH (to correctly deprotonate or protonated the protein), and then add maybe the ion of interest (in our case Gd3+) before neutralising the enironment (the overall charge of the system need to be 0) by adding counterions (Na^+ and Cl^-).
> One of the important steps of the system preparation is the choice of the force field and of the water model. Those are theoretical model set up to describe the interation between the different particles in our system. Choosing the correct system is determinant for the simulation to be correct.


### The 12-6-4 LJ non bondd model
The 12-6-4 LJ non-bonded model allows for a better representation of the interaction between highly charged metal ions and proteins, consequently it is important to use it when studying protein complex.
#### File preparation
To start using amber mit metal complex and the 12-6-4 model, the following data and files are needed. 
- The *.pdb* file of the protein.
- The location of the metal ion in the protein.
- The different input files:
  - min_1.in
  - min_2.in
  - heat.in
  - ...
- The leap.in file (it has some specificity compared to a normal leap.in file)
Then we need to generate **two** important files: the *mol2* and *frcmod* files. Fo this we use **antechamber** and **parmchk**.
First we generate the mol2 file with the following command:
```
#if amber is not already loaded
module load Amber
antechamber -i MOL.pdb -f pdb -o MOL.mol2 -fo mol2 -c bcc -nc 0 -m 1
```
The *-nc* variable is the charge of the molecule and the *-m* is the multiplicity.
Then when the mol2 file done we can use it to form the frcmod file:
```
parmchk - MOL.mol2 -f mol2 -o MOL.frcmod -fo 
```
Now we have all what we need,let's run **tleap**:
```
tleap -s -f leap.in
```
> The tleap.in file needs to possess the addC4Type line to do have the C4 term added.
