# GROMACS Tutorial
This page is dedicated to learning the basis of GROMACS. GROMACS is programm calculating different parameter of molecules in a given environement in order to get the final calculation. To "_feed_" the program multiple saller fiel are needed. Each one of those files are filled with the information about either your molecule or it environement. Consequently to change a parameter of the simulation or to share the simulation only those small files are necessary.

## The structure of GROMACS
GROMACS has a fixed strutured, an order in which different operations are done. At each step, something new is calculated. Schematically the process is devided into the following __seven__ steps/operations:
1. __The generation of the topology__
2. __The defenition of the box and solvate__
3. __The ions addition__
4. __The energy minimization__
5. __The equilibration__
6. __The production of the MD (_Molecular Dynamic_)__
7. __The analysis__

### The generation of the topology
#### What is the _topology_
The molecular topology is the description of the atom _and_ the combination of atom __on wich_ the potential functions are going to be applied on. In GROMACS the topology file is the __.top__ extension. It formed at the start, because every other calculation will need to be done based on it (For further detail on the interaction function adn the force field, see force field). Be aware that GROMACS is biological system oriented sofware, meaning that only the atoms present in biological system are parameterized in the force field. It is important to note that the __velocity__ and the __trajectory__ of the atoms are not stored in the .top file but rather in the __.gro__ (position and velocites) and __.trr__(positions, velocities and force). 
### The defenition of the box and solvate
### The ions addition
### The energy minimization
### The equilibration
### The production of the MD 
### The analysis
