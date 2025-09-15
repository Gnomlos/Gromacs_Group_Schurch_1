#!/usr/bin/env python3
# This script is used to map the protein with their superfamily of CATH3D id (funfam).
# Also please make sure that the libraries are installed 

import argparse, time, requests, pandas as pd
from pathlib import Path

PDBe = "https://www.ebi.ac.uk/pdbe/api/mappings"
UniProt = "https://rest.uniprot.org/uniprotkb"

def get_json(url, retries=3):
    for _ in range(retries):
        r = requests.get(url, headers={"User-Agent":"gd3-cath-mapper"})
        if r.status_code == 200:
            return r.json()
        time.sleep(1)
    return None

def pdb_to_uniprot(pdb_id):
    """Return list of (uniprot_id, protein_name) for a PDB."""
    data = get_json(f"{PDBe}/{pdb_id.lower()}")
    if not data or pdb_id.lower() not in data: return []
    out = []
    for acc, info in data[pdb_id.lower()].get("UniProt",{}).items():
        name = info.get("name") or info.get("entry_name") or ""
        out.append((acc, name))
    return out

def uniprot_to_cath(acc):
    """Return list of (cath_id) from UniProt entry."""
    data = get_json(f"{UniProt}/{acc}.json")
    if not data: return []
    cath = []
    for xref in data.get("uniProtKBCrossReferences",[]):
        if xref.get("database")=="Gene3D":
            props = {p["key"]:p["value"] for p in xref.get("properties",[])}
            cath_id = props.get("CATH id") or xref.get("id")
            #cath_name = props.get("Entry name") or props.get("Entry description")
            if cath_id:
                cath.append((cath_id))
    return cath

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--pdb-list",required=True,help="CSV with pdb_id column")
    ap.add_argument("--input-table",help="Original gd3_ranking.csv to annotate")
    ap.add_argument("--protein-col",default="protein",help="Column in input-table with pdb-coded protein names")
    ap.add_argument("--out-prefix",default="gd3_cath")
    args=ap.parse_args()

    pdb_ids=pd.read_csv(args.pdb_list)["pdb_id"].dropna().astype(str).str.lower().unique()

    rows=[]
    for i,pdb in enumerate(pdb_ids,1):
        for acc,name in pdb_to_uniprot(pdb):
            for cid in uniprot_to_cath(acc):
                rows.append({"pdb_id":pdb,"uniprot_id":acc,"protein_name":name,
                             "cath_superfamily_id":cid})
        if i%50==0: print(f"...{i}/{len(pdb_ids)} done")
        time.sleep(0.1)

    cath_df=pd.DataFrame(rows)
    cath_df.to_csv(f"{args.out_prefix}_per_protein.csv",index=False)
    if args.input_table:
        tbl = pd.read_csv(args.input_table)

        # Extract the 4-char PDB ID from the protein column (handles "3B2D_protonated_out")
        tbl["_pdb_id"] = (
            tbl[args.protein_col]
            .astype(str)
            .str.extract(r"^([0-9A-Za-z]{4})", expand=False)
            .str.lower()
    )

    # Merge on pdb_id
    annotated = tbl.merge(cath_df, left_on="_pdb_id", right_on="pdb_id", how="left")

    out_tbl = Path(f"{args.out_prefix}_annotated.csv")
    annotated.to_csv(out_tbl, index=False)
    print(f"Wrote annotated table to {out_tbl}")
    # --- Summary statistics ---
    n_proteins = cath_df["uniprot_id"].nunique(dropna=True)
    n_superfam = cath_df["cath_superfamily_id"].nunique(dropna=True)

    summary_text = (
        f"There are {n_proteins} unique proteins "
        f"repartitioned over {n_superfam} CATH superfamilies."
    )
    print(summary_text)

    # Save to a small text file
    summary_path = Path(f"{args.out_prefix}_summary.txt")
    with open(summary_path, "w") as f:
        f.write(summary_text + "\n")
    print(f"Wrote summary to {summary_path}")

if __name__=="__main__":
    main()
