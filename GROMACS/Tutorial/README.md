# GROMACS Tutorial

This page is dedicated to learning the basis of GROMACS. GROMACS is programm calculating different parameter of molecules in a given environement in order to get the final calculation. To "_feed_" the program multiple saller fiel are needed. Each one of those files are filled with the information about either your molecule or it environement. Consequently to change a parameter of the simulation or to share the simulation only those small files are necessary.

## The structure of GROMACS

GROMACS has a fixed strutured, an order in which different operations are done. At each step, something new is calculated. Schematically the process is devided into the following **seven** steps/operations:

1.  **The generation of the topology**
    
2.  **The definition of the box and solvate**
    
3.  **The ions addition**
    
4.  **The energy minimization**
    
5.  **The equilibration**
    
6.  **The production of the MD (_Molecular Dynamic_)**
    
7.  **The analysis**
    

### The generation of the topology

#### What is the _topology_

The molecular topology is the description of the atom _and_ the combination of atom \__on wich_ the potential functions are going to be applied on. In GROMACS the topology file is the **.top** extension. It formed at the start, because every other calculation will need to be done based on it (For further detail on the interaction function adn the force field, see force field). Be aware that GROMACS is biological system oriented sofware, meaning that only the atoms present in biological system are parameterized in the force field. It is important to note that the **velocites** and the **trajectories** of the atoms are not stored in the .top file but rather in the **.gro** (position and velocites) and **.trr**(positions, velocities and force).

GROMACS recognized **three** particle types:

| Particle | Symbole |
| --- | --- |
| atom | A |
| shell | S |
| virtual site | V (or D) |

Shells are only used for the polarizable models.

##### Atom types

Each atom is defined into a foce field, by its mass and name. This "_list_" is found in a file called the _atomtypes.atp_. (.atp) If a **new** atom type need to be added, it is in this file that it needs to be done. Be aware that only the command gmx pdb2gmx use this file. The interactions described in the .atp and \[ atomtypes\] section of the topology file describe only **non bonding interaction** like the Van der Walls interactions.

##### Virtual site

Virtual sites are interactio site constructed based on the position of other particle positions on which certain interactions are located.

Virtual site can be set up in the .top file under the \[ virtual\_sites? \] section with the **?** standing for the number of virtual site present. This being a rather complicated subject, it is recommander to read the section dedicated to it on [https://manual.gromacs.org/current/reference-manual/topologies/particle-type.html](https://manual.gromacs.org/current/reference-manual/topologies/particle-type.html) .

#### Parameter files

Atom have different proprieties which are called, _parameters_.

| Propriety | Symbol | Unit |
| --- | --- | --- |
| Type (Name of the atom) | \- | \- | 
| Mass | m | a.m.u |  
| Charge | q | electron |  
| Epsilon | \epsilon | kJ/mol |
| sigma | \sigma | nm |

Those proprieties are initialy stored in the `atomtypes.atp` for the atom mass, the charge are found in the __.rtp__ file, where the \epsilon (denoted as W or c12)and \sigma (denoted as V or c6) paramters are found in the `ffnonbonded.itp`. The different proprieties are then combine in the topology file. 
> It is to be noted that the \sigma and \epsilon factors are important in the calculation of the __non-bonded__ interaction, using the [Lennard-Jone interaction](https://manual.gromacs.org/current/reference-manual/functions/nonbonded-interactions.html#lj).

Bonded parameter are found in a separate file:`ffbonded.itp` which are then inclueded in the topology file under either the `[ bonds ]`,  `[ angle ]`,`[ dihedral ]`,... section.
#### Molecule definition
Until now, only the ocntext of the atom was discussed, __but__ GROMACS is used to simulate whole molecules. In the topology file, the `[ moleculetype]` entry allows to structure in term of molecules, allowing the have only one single entry for every moleucle type present in our simulation --> this implies that in the topology file all the molecule present in the reaction are present. This represent alo a challenge since GROMACS needs to be told how those molecule are interacting with each other, in other term how does the same moleucle behave when it meets the itself. This is done in the _mdp_ file.
> Molecule  which convalently bond a ligand need to be describe in the `[ moleculetype]` entry. If not convalently bonded then they both needs their own separate entry in `[ moleculetype]`. The non-convalent interaction are define at the end of the topology file in the `[ intermolecular_interactions ]` entry.
> Another entry of the topology file is the entry `[ pair ]` which allows for the specific definition of electrostatic interaction between a pair of atom. It is to be noted that interactions do not take plpace between atom which are to close to each other in the molecule. This is the concept of __exclusion__. Extra exclusion can be added in the `[ exclusion ]` entry at the end of the topology file.



### The definition of the box and solvate

### The ions addition

### The energy minimization

### The equilibration

### The production of the MD

### The analysis
