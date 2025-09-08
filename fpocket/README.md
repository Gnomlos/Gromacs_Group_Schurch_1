# Fpocket Toolkit

The **fpocket** suite is an open-source set of programs designed to detect and analyze binding pockets in proteins. It is widely used in structural bioinformatics and computational drug discovery to identify potential ligand-binding sites and to evaluate their physicochemical properties.  

---

## Table of Contents

- [Installation](#installation)
- [Components](#components)
- [Basic Usage](#basic-usage)
- [Fpocket Suite Parameters](#fpocket-suite-parameters)
  - [fpocket](#fpocket)
  - [dpocket](#dpocket)
  - [mdpocket](#mdpocket)
  - [tpocket](#tpocket)
- [Fpocket Output Overview](#fpocket-output-overview)
- [Notes](#notes)
- [Citation](#citation)

---

## Installation

Clone the repository and build the executables:  

```bash
git clone https://github.com/Discngine/fpocket.git
cd fpocket
make
```

This will generate the main executables (`fpocket`, `mdpocket`, `dpocket`, and `tpocket`) in the `bin/` directory.  
To run them easily, add the binaries to your `PATH`:  

```bash
export PATH=$PATH:/path/to/fpocket/bin
```

---

## Components

- **`fpocket`**  
  Detects potential binding pockets in protein structures (PDB files) using a Voronoi tessellation and alpha-sphere method.

- **`mdpocket`**  
  Analyzes molecular dynamics trajectories to study the stability and evolution of pockets over time.

- **`dpocket`**  
  Scores and characterizes pockets, focusing on druggability and detailed properties. Typically used after `fpocket`.

- **`tpocket`**  
  Provides statistical analysis and clustering of pockets, useful for comparisons across proteins or conditions.

---

## Basic Usage

Run **fpocket** on a PDB structure:

```bash
fpocket -f protein.pdb
```

- `-f protein.pdb` : Input protein structure (PDB format).  

Output: a folder named `<protein>_out/` containing:  
- `pockets/` : individual PDB files for each pocket  
- `descriptors.txt` : pocket descriptors and scoring values  
- `out.pdb` : protein structure annotated with detected pockets  

---

# Fpocket Suite Parameters

Below is an overview of the most important parameters for each tool in the fpocket suite.  
For full details, always check the built-in help (`<tool> -h`).  

---

## `fpocket`

Detects potential binding pockets in a static protein structure.

| Parameter | Description | Default |
|-----------|-------------|---------|
| `-f <file>` | Input protein structure in **PDB** format (required). | – |
| `-p <path>` | Output directory prefix. Results written to `<prefix>_out/`. | Based on input filename |
| `-s <min>` | Minimum number of alpha spheres per pocket. | 50 |
| `-m <max>` | Maximum number of alpha spheres per pocket. | 1000 |
| `-d <dist>` | Maximum distance (Å) between alpha spheres for clustering. | 4.5 |
| `-r <radius>` | Minimum alpha sphere radius. | 3.0 |
| `-c` | Enable clustering of nearby pockets. | Off |
| `-t` | Output pocket PDB files with transparency attributes. | Off |
| `-h` | Show help message. | – |

**Example:**

```bash
fpocket -f protein.pdb -s 100 -d 3.5
```

---

## `dpocket`

Scores and characterizes pockets detected by `fpocket`.

| Parameter | Description | Default |
|-----------|-------------|---------|
| `-f <file>` | Input protein structure (PDB). | – |
| `-p <path>` | Path to fpocket output directory. | – |
| `-i <id>` | Analyze only the pocket with given ID. | – |
| `-h` | Show help message. | – |

**Example:**

```bash
dpocket -f protein.pdb -p protein_out/
```

---

## `mdpocket`

Analyzes molecular dynamics (MD) trajectories to monitor pocket evolution.

| Parameter | Description | Default |
|-----------|-------------|---------|
| `-f <topology>` | Reference protein topology (PDB or PQR). | – |
| `-d <traj>` | MD trajectory file (DCD, TRR, etc.). | – |
| `-o <out>` | Output directory prefix. | `mdpocket_out/` |
| `-s <stride>` | Stride (frames to skip) for trajectory analysis. | 1 |
| `-t <time>` | Maximum time to analyze (in ps). | Full trajectory |
| `-h` | Show help message. | – |

**Example:**

```bash
mdpocket -f protein.pdb -d traj.dcd -s 10 -o md_analysis
```

---

## `tpocket`

Performs statistical analysis and clustering across multiple pocket datasets.

| Parameter | Description | Default |
|-----------|-------------|---------|
| `-l <list>` | List of fpocket output directories to compare. | – |
| `-o <out>` | Output file prefix. | `tpocket_out` |
| `-c` | Cluster pockets across datasets. | Off |
| `-h` | Show help message. | – |

**Example:**

```bash
tpocket -l pocket_list.txt -o comparison
```

---

## Fpocket Output Overview

When running `fpocket`, the program generates an output folder named `<protein>_out/`.  
This contains:

- **`out.pdb`** – input protein annotated with pocket spheres.  
- **`pockets/`** – directory containing PDB files for each detected pocket.  
- **`descriptors.txt`** – numerical descriptors (volume, polarity, hydrophobicity, druggability score, etc.) for each pocket.  
- **`ranked_pockets.pdb`** – pockets sorted by fpocket scoring function.  
- **`info.txt`** – summary of detected pockets and scoring details.  

These files allow both visualization (in PyMOL, Chimera, etc.) and quantitative analysis of potential binding sites.

---

## Notes

- All tools support `-h` for help and usage.  
- File formats must follow the PDB/MD conventions supported by the suite.  
- Typical workflow:  
  1. Detect pockets with `fpocket`  
  2. Characterize with `dpocket`  
  3. Monitor over MD with `mdpocket`  
  4. Compare across runs with `tpocket`  

---

## Citation

If you use **fpocket** in your research, please cite:  

> Le Guilloux, V., Schmidtke, P., & Tuffery, P. (2009). Fpocket: An open source platform for ligand pocket detection. *BMC Bioinformatics, 10*(1), 168. https://doi.org/10.1186/1471-2105-10-168  

### BibTeX

```bibtex
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
```
