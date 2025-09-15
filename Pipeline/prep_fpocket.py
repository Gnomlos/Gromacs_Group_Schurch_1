#!/usr/bin/env python3
# It is important to note that this script imply that fpocket is installed and that all the libraries are also intalled.
# So this script will do the folowing:
#  1. Clean the pdb file of any water or other atom/molecules whcih does not belong to the protein.
#  2. Protonate/deprotone the protein accordingly to the pH of interest, in our case pH = 7.4
#  3. Run fpocket automatically and save the results, with a folder by proteins

import argparse
import pathlib
import subprocess
from tqdm import tqdm
from pdbfixer import PDBFixer
from openmm.app import PDBFile

def clean_pdb(input_pdb: pathlib.Path, out_dir: pathlib.Path) -> pathlib.Path:
    """Remove non-ATOM/HETATM records and add END."""
    out_path = out_dir / input_pdb.name
    with open(input_pdb, "r") as f_in, open(out_path, "w") as f_out:
        for line in f_in:
            if line.startswith(("ATOM", "HETATM")):  # keep ligands/ions
                f_out.write(line)
        f_out.write("END\n")
    return out_path


def protonate_pdb(input_pdb: pathlib.Path, out_dir: pathlib.Path, ph: float) -> pathlib.Path:
    """Add missing hydrogens at a given pH using PDBFixer."""
    fixer = PDBFixer(filename=str(input_pdb))

    # Must be called in order
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(ph)

    out_path = out_dir / (input_pdb.stem + "_protonated.pdb")
    with open(out_path, "w") as f_out:
        PDBFile.writeFile(fixer.topology, fixer.positions, f_out, keepIds=True)

    return out_path


def run_fpocket(pdb_path: pathlib.Path, out_dir: pathlib.Path):
    """Run fpocket on a PDB file and collect status."""
    subprocess.run(
        ["fpocket", "-f", str(pdb_path), "-m", "4.6", "-M", "7.5"],
        cwd=out_dir,
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL
    )
    return pdb_path.name, "OK"


def process_pipeline(pdb_file: pathlib.Path, args, log_file: pathlib.Path):
    """Full pipeline: clean → protonate → run fpocket, with error handling."""
    cleaned_dir = args.output / "cleaned"
    prot_dir = args.output / "protonated"
    fpocket_dir = args.output / "fpocket_results"

    for d in (cleaned_dir, prot_dir, fpocket_dir):
        d.mkdir(parents=True, exist_ok=True)

    try:
        # Step 1: Clean
        cleaned_pdb = clean_pdb(pdb_file, cleaned_dir)

        # Step 2: Protonate
        prot_pdb = protonate_pdb(cleaned_pdb, prot_dir, args.ph)

        # Step 3: Run fpocket
        name, status = run_fpocket(prot_pdb, fpocket_dir)

    except Exception as e:
        name, status = pdb_file.name, f"ERR:{type(e).__name__} {e}"
        with open(log_file, "a") as log:
            log.write(f"{pdb_file.name}: {status}\n")

    return name, status


def main():
    parser = argparse.ArgumentParser(description="Batch clean, protonate, and run fpocket (sequential, robust)")
    parser.add_argument("input", help="Input PDB file or folder")
    parser.add_argument("-o", "--output", type=pathlib.Path, default=pathlib.Path("results"),
                        help="Output folder (default: ./results)")
    parser.add_argument("--ph", type=float, default=7.4, help="Target pH (default: 7.4)")
    args = parser.parse_args()

    input_path = pathlib.Path(args.input)
    if input_path.is_dir():
        pdb_files = [p for p in sorted(input_path.glob("*.pdb")) if "_out" not in p.name]
    else:
        pdb_files = [input_path]

    args.output.mkdir(parents=True, exist_ok=True)
    log_file = args.output / "error_log.txt"

    results = []
    for f in tqdm(pdb_files, desc="Processing", unit="file"):
        results.append(process_pipeline(f, args, log_file))

    # Print results summary
    success = sum(1 for _, status in results if status == "OK")
    failed = len(results) - success

    print("\n=== SUMMARY ===")
    print(f"✅ Successful: {success}")
    print(f"❌ Failed:     {failed}")
    if failed > 0:
        print(f"See error log: {log_file}")


if __name__ == "__main__":
    main()
