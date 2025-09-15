# run_fpocket_batch.py
# So this scrip is mainly use for runing f-pocket alone over a range of protein. Beware of the following: it clean the data by removing possible non-protein atom, like water. But it DOESN'T adapt the protein to the pH. 

import subprocess, os, pathlib
from tqdm import tqdm   # For the progress bar

PDB_DIR = pathlib.Path("/home/gadolinium/Documents/Protein_pipeline/pdb_structures/") # Here change the path where the pdb structure are downloaded.
PDBS = [p for p in sorted(PDB_DIR.glob("*.pdb")) if "_out" not in p.name]

def cleaning(path_pdb, o_dir):
    cleaned = os.path.join(o_dir, "cleaned") # If you want to change the export name, replace "cleaned" with what you want
    os.makedirs(cleaned, exist_ok=True)
    cl_pbs_list = []

    # tqdm loop for cleaning
    for pdb in tqdm(path_pdb, desc="Cleaning PDBs", unit="file"):
        out_path = os.path.join(cleaned, pathlib.Path(pdb).name)
        with open(pdb, "r") as f_in, open(out_path, "w") as f_out:
            for line in f_in:
                if line.startswith(("ATOM", "HETATM")): #sometime it is good to supress "HETATM" because it can pause some issue latter on.
                    f_out.write(line)
            f_out.write("END\n")  # ensure END is present
        cl_pbs_list.append(pathlib.Path(out_path))

    return cl_pbs_list

def run_fpocket(pdb_path: pathlib.Path):
    try:
        subprocess.run(
            ["fpocket", "-f", str(pdb_path), "-m", "4.6", "-M", "7.5"], # Here change the fpocket parameter you wish to use.
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL
        )
        return pdb_path.name, "OK"
    except subprocess.CalledProcessError as e:
        return pdb_path.name, f"ERR:{e.returncode}"

if __name__ == "__main__":
    # Cleaning phase with tqdm
    clean_pbs = cleaning(PDBS, PDB_DIR) 
    
    # Fpocket phase with tqdm
    for pdb in tqdm(clean_pbs, desc="Running fpocket", unit="file"):
        name, status = run_fpocket(pdb)
        tqdm.write(f"{name}: {status}")
