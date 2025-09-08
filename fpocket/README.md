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
### Fpocket Parameters

The `fpocket` program accepts several parameters to control pocket detection.  
Here are the most commonly used options:

| Parameter | Description | Default |
|-----------|-------------|---------|
| `-f <file>` | Input protein structure in **PDB** format (required). | – |
| `-p <path>` | Output directory prefix. Results will be written to `<prefix>_out/`. | Based on input filename |
| `-s <min>` | Minimum number of alpha spheres per pocket. Filters out very small pockets. | 50 |
| `-m <max>` | Maximum number of alpha spheres per pocket. | 1000 |
| `-i <id>` | Keep only the pocket with the specified ID (useful for re-analysis). | – |
| `-d <dist>` | Maximum distance (in Å) between alpha spheres to be clustered into the same pocket. | 4.5 |
| `-r <radius>` | Minimum radius of alpha spheres considered. | 3.0 |
| `-c` | Enable clustering of close pockets into larger ones. | Off |
| `-t` | Output pocket PDB files with transparency attributes (for visualization). | Off |
| `-h` | Display the help message with all available options. | – |

---

### Example: Custom Parameters

Run fpocket with stricter filtering (only larger pockets with at least 100 alpha spheres):

```bash
fpocket -f protein.pdb -s 100
```
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
