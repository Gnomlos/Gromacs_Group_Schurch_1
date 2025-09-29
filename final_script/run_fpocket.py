#!/usr/bin/env python3
import argparse
import pathlib
import subprocess
from typing import Dict, List, Tuple

import pandas as pd
from tqdm import tqdm

from pdbfixer import PDBFixer
from openmm.app import PDBFile
from Bio.PDB import PDBParser, PPBuilder
from Bio.SeqUtils import IsoelectricPoint as IP
from Bio.SeqUtils import ProtParam

# run fpocket but also clean the pdb of possible cofactor and also protonate the pdb structure accordingly to the wanted pH (genaraly pH=7.4)

# --------- PDB processing ----------

def clean_pdb(input_pdb: pathlib.Path, out_dir: pathlib.Path) -> pathlib.Path:
    """
    Keep only protein ATOM records (drop all HETATM, incl. waters/ligands/ions).
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / input_pdb.name
    with open(input_pdb, "r") as f_in, open(out_path, "w") as f_out:
        for line in f_in:
            # canonical PDB: protein atoms are "ATOM"
            if line.startswith("ATOM"):
                f_out.write(line)
        f_out.write("END\n")
    return out_path


def protonate_pdb(input_pdb: pathlib.Path, out_dir: pathlib.Path, ph: float) -> pathlib.Path:
    """
    Add missing atoms and hydrogens at target pH using PDBFixer.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    fixer = PDBFixer(filename=str(input_pdb))
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(ph)

    out_path = out_dir / (input_pdb.stem + "_prot.pdb")
    with open(out_path, "w") as f_out:
        PDBFile.writeFile(fixer.topology, fixer.positions, f_out, keepIds=True)
    return out_path


def run_fpocket(pdb_path: pathlib.Path, out_dir: pathlib.Path) -> Tuple[str, str]:
    """
    Run fpocket on a PDB file. Returns (basename, status).
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        ["fpocket", "-f", str(pdb_path)],
        cwd=out_dir,
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL
    )
    return pdb_path.name, "OK"


# --------- Charge computation from PDB sequence ----------

def sequences_from_pdb(pdb_file: pathlib.Path) -> Dict[str, str]:
    """
    Parse sequences per chain from a PDB using Biopython.
    Returns dict: chain_id -> sequence (one-letter codes).
    Notes:
      - Uses ATOM records (after cleaning), so only protein chains remain.
      - Missing residues in coordinates are absent from the sequence.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("prot", str(pdb_file))
    ppb = PPBuilder()
    chain_seqs: Dict[str, str] = {}

    # Build peptides for each chain; PPBuilder returns continuous fragments.
    for model in structure:
        for chain in model:
            frags = ppb.build_peptides(chain)
            if not frags:
                continue
            seq = "".join(str(pp.get_sequence()) for pp in frags)
            if seq:
                chain_id = chain.id if chain.id is not None else "?"
                chain_seqs[chain_id] = seq
        break  # use first model if multiple
    return chain_seqs


def charge_of_sequence(seq: str, ph: float = 7.0) -> float:
    """
    Compute formal charge of a peptide sequence at a given pH using IsoelectricPoint.
    """
    aa_counts = ProtParam.ProteinAnalysis(seq).count_amino_acids()
    # IsoelectricPoint expects the raw sequence and an aa count dict
    ip = IP.IsoelectricPoint(seq, aa_counts)
    return ip.charge_at_pH(ph)


def compute_charge_from_pdb(pdb_file: pathlib.Path, ph: float = 7.0) -> Tuple[float, Dict[str, float]]:
    """
    Compute total and per-chain formal charge from a (protonated or plain) PDB.
    Returns (total_charge, {chain_id: charge}).
    """
    chains = sequences_from_pdb(pdb_file)
    per_chain: Dict[str, float] = {}
    for cid, seq in chains.items():
        if seq:  # skip empty
            try:
                per_chain[cid] = charge_of_sequence(seq, ph)
            except Exception:
                per_chain[cid] = float("nan")
        else:
            per_chain[cid] = float("nan")
    total = sum(c for c in per_chain.values() if pd.notna(c)) if per_chain else float("nan")
    return total, per_chain


# --------- Orchestration ----------

def main():
    parser = argparse.ArgumentParser(
        description="Clean → protonate → fpocket and append formal charge (from PDB sequence) to mapping CSV."
    )
    parser.add_argument("mapping_csv", help="Input CSV (e.g., uniprot_pdb_cath_mapping.csv)")
    parser.add_argument("-i", "--input", type=pathlib.Path, required=True,
                        help="Folder containing input *.pdb files")
    parser.add_argument("-o", "--output", type=pathlib.Path, default=pathlib.Path("prep_results"),
                        help="Output folder (default: ./prep_results)")
    parser.add_argument("--ph", type=float, default=7.0, help="Target pH for protonation & charge (default: 7.0)")
    args = parser.parse_args()

    df = pd.read_csv(args.mapping_csv)
    if "PDB_ID" not in df.columns:
        raise ValueError("mapping_csv must contain a 'PDB_ID' column.")

    clean_dir = args.output / "cleaned"
    prot_dir = args.output / "protonated"
    fpocket_dir = args.output / "fpocket_results"
    args.output.mkdir(parents=True, exist_ok=True)

    # Cache charges per PDB to avoid recomputation when multiple rows share the same PDB_ID
    charge_cache: Dict[str, Tuple[float, Dict[str, float]]] = {}

    total_charges: List[float] = []
    chain_charge_strs: List[str] = []

    for _, row in tqdm(df.iterrows(), total=len(df), desc="Processing", unit="row"):
        pdb_id = str(row["PDB_ID"]).strip()
        pdb_in = args.input / f"{pdb_id}.pdb"

        if not pdb_in.exists():
            total_charges.append(float("nan"))
            chain_charge_strs.append("")
            continue

        if pdb_id not in charge_cache:
            try:
                # 1) Clean (drop all HETATM / waters / ligands)
                cleaned = clean_pdb(pdb_in, clean_dir)

                # 2) Protonate at requested pH
                prot = protonate_pdb(cleaned, prot_dir, args.ph)

                # 3) Run fpocket on protonated structure
                run_fpocket(prot, fpocket_dir)

                # 4) Compute charge from the (protonated) PDB sequence(s)
                total_c, per_chain = compute_charge_from_pdb(prot, ph=args.ph)

                # store
                charge_cache[pdb_id] = (total_c, per_chain)

            except Exception as e:
                charge_cache[pdb_id] = (float("nan"), {})

        total_c, per_chain = charge_cache[pdb_id]
        total_charges.append(total_c)
        # Store per-chain charges as a compact string like "A:1.2;B:-3.8"
        if per_chain:
            chain_str = ";".join(f"{cid}:{per_chain[cid]:.4g}" for cid in sorted(per_chain))
        else:
            chain_str = ""
        chain_charge_strs.append(chain_str)

    # Append to dataframe
    df["Total_Charge_pH7"] = total_charges if abs(args.ph - 7.0) < 1e-6 else pd.NA
    df[f"Total_Charge_pH{args.ph:g}"] = total_charges
    df[f"Chain_Charges_pH{args.ph:g}"] = chain_charge_strs

    out_csv = args.output / "mapping_with_charge_from_pdb.csv"
    df.to_csv(out_csv, index=False)
    print(f"\nUpdated mapping saved to: {out_csv}")


if __name__ == "__main__":
    main()
