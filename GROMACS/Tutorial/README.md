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

| Propriety | Symbol | Unit |  |
| --- | --- | --- | --- |
| Type (Name of the atom) | \- | \- |  |
| Mass | Cell 2-2 | Cell 2-3 |  |
| Charge | Cell 3-2 | Cell 3-3 |  |
| Epsilon | Cell 4-1 | Cell 4-2 | Cell 4-3 |
| sigma | Cell 5-2 | Cell 5-3 |  |

In the `[ atomtypes ]` section of your .top file it

##### The non bonded interactions

Non bonded interactions represent the interactions between non-bonded atoms. They are represented by the following equations: $$

### The definition of the box and solvate

### The ions addition

### The energy minimization

### The equilibration

### The production of the MD

### The analysis