# ORCA -- An introduction to quantum chemistry computing

In the context of being able to "*predict*" the possible bonding point of $Gd^{3+}$ and $Gd-DOTA$ to protein using a high throughput method, it important to be able to calculate the affinity of the protein with the gadolinium. To find this affinity, the 12-6-4 potentials between the gadolinium and the different amino acid needs to calculated. 

> # 12-6-4 potential
> The 12-6-4 potentials are **parameters** used in the following equation:
> $aa$
> We need those parameter for every possible interactions.

## Single point calculation

Single point calculations calculate the total energy of a system at a given geometry. Since those calculation are complicated and expensiv to manually do, a quantum chemiststry software is needed. In our lab the software ORCA was used due to its simplicity of installation and use. The goal of the those calculation, as previously stated, is to find the 12-6-4 potentials of Gd and amino-acids(AA). In order to achieve this we need to calculate the **interaction energy** between the metal and the amino acid at *different distances*. To isolate the interaction energy of the Gd-AA from the energy of the structure we need to calculate the enrgy of system of the Gd and the AAs alone, then of both of them at different distances. The following equation allows us then to calculate the interaction energy:<br>
$E_{Gd-AA} = E_{sys} - E_{Gd} - E_{AA}$<br>

