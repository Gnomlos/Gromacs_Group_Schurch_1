# FPOCKET
The fpocket suite is an open-source set of programs designed to detect and analyze binding pockets in proteins. It is to identify potential ligand-binding sites and to evaluate their physicochemical properties.
## Installation
So fpocket is only available on linux (if I am not mistaken). To install it, clone the repository and build the executables:

```bash
git clone https://github.com/Discngine/fpocket.git
cd fpocket
make
```
Then update the path:
```bash
export PATH=$PATH:/path/to/fpocket/bin
```
> Note: It is import to note that you need to be admin to install the program.

## Tools
Fpocket is actually 4 different tools grouped in one package: fpocket, mdpocket, dpocket and tpocket. Each one have its own utility. In the context of our work we only use (until now) fpocket. However I will briefly described the role of each one.
### fpocket
### mdpocket
### dpocket
### tpocket

## Syntax

## Reference
@article{LeGuilloux2009Fpocket,
  title   = {Fpocket: An open source platform for ligand pocket detection},
  author  = {Le Guilloux, Vincent and Schmidtke, Peter and Tuffery, Pierre},
  journal = {BMC Bioinformatics},
  volume  = {10},
  number  = {1},
  pages   = {168},
  year    = {2009},
  doi     = {10.1186/1471-2105-10-168}
}
