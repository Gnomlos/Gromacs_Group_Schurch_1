# ORCA -- An introduction to quantum chemistry computing

In the context of being able to "*predict*" the possible bonding point of $Gd^{3+}$ and $Gd-DOTA$ to protein using a high throughput method, it important to be able to calculate the affinity of the protein with the gadolinium. To find this affinity, the 12-6-4 potentials between the gadolinium and the different amino acid needs to calculated. 

> # 12-6-4 potential
> The 12-6-4 potentials are **parameters** used in the following equation:
> $aa$
> We need those parameter for every possible interactions.

## Single point calculation

Single point calculations calculate the total energy of a system at a given geometry. Since those calculation are complicated and expensiv to manually do, a quantum chemiststry software is needed. In our lab the software ORCA was used due to its simplicity of installation and use. The goal of the those calculation, as previously stated, is to find the 12-6-4 potentials of Gd and amino-acids(AA). In order to achieve this we need to calculate the **interaction energy** between the metal and the amino acid at *different distances*. To isolate the interaction energy of the Gd-AA from the energy of the structure we need to calculate the enrgy of system of the Gd and the AAs alone, then of both of them at different distances. The following equation allows us then to calculate the interaction energy:<br>
$E_{Gd-AA} = E_{sys} - E_{Gd} - E_{AA}$<br>
It is to be noted that the structure of the AA needs to be optimized before the energy is calculated, see following section. In orca the following code is given to calculate the single point energy:
```
!HF DEF2-SVP
* xyz 0 1
O         -3.56626        1.77639        0.00000
H         -2.59626        1.77639        0.00000
H         -3.88959        1.36040       -0.81444
*
```
The first line gives the method used to calculation (see the method section for more information on which method to use). The `!` is obligatory no matter what method is used. The next section of the code represent, the studied molecule. It is basically the corrdinate of the different atom, using the classical syntax of a `.xyz` file. It is important to note that in the case of a complex, the Gd-AA complex for example. All coordinate needs to be put in a ***single section*** as the the software calculate the energy of a system. The 2 number after the `xyz` in the code are the charge and the multiplicity of the system. This need to be specified for every system. 




## Geometry optimization

Molecules geometries have intresic energy which, in most case, are at the lower in the conformation at wich the molecule is found. So to calculate the energy accuratly the geometry need to be "*optimized*", meaning at it lowest energy. 

