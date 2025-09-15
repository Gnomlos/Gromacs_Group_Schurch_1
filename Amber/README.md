# AMBER: A Molecular dynamic simulation software
Amber, much like GROMACS is one of the most used md software. The main difference with GROMACS is however its **flexibility**. Over time multiple settings and tools for different case were devlopped. In our case, their is two main tools
wich are very intersesting: MCBP.py and MMPBSA.py. Both of wich will serve us to simulate the Gd ions and to calculate the bonding energy of this ion to the protein (possibly molecule). 
> For the installation of Amber, please consult the instruction in the *Amber_source_file* folder.
## A quick introduction into molecular dynamic
Before starting using Amber it is preferable to understand what exactly we are doing. This will helps us to better interpret our results and check for any problem as well as correcting problems. **Molecular dynamic simulation** is *simply*
the simulation of a chemical environement using a computer. However contrary to quantum chemistry computin (like **ORCA**), it is based on classical physics equation and not on quantum theory. The advantage of this is that it allows for a 
simpler, lighter software which gives overall good results; it is however a bit less precise. 

Ok, but how do we simulate our experience? Simply buy telling the software what are the conditions of our experiment **in a precise order**. The order of steps is important and it is actually where the major errors can come from; so lets decompose the different steps and explain what each do. 
