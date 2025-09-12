#!/usr/bin/env python3
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
        name = info.get("protein_name") or info.get("entry_name") or ""
        print(acc," ",name)
        out.append((acc, name))
    return out

def uniprot_to_cath(acc):
    """Return list of (cath_id, cath_name) from UniProt entry."""
    data = get_json(f"{UniProt}/{acc}.json")
    if not data: return []
    cath = []
    for xref in data.get("uniProtKBCrossReferences",[]):
        if xref.get("database")=="Gene3D":
            props = {p["key"]:p["value"] for p in xref.get("properties",[])}
            cath_id = props.get("CATH id") or xref.get("id")
            cath_name = props.get("Entry name") or props.get("Entry description")
            if cath_id:
                cath.append((cath_id,cath_name))
                print("cath ", cath_id)
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
            for cid,cname in uniprot_to_cath(acc):
                rows.append({"pdb_id":pdb,"uniprot_id":acc,"protein_name":name,
                             "cath_superfamily_id":cid,"cath_superfamily_name":cname})
        if i%50==0: print(f"...{i}/{len(pdb_ids)} done")
        time.sleep(0.1)

    cath_df=pd.DataFrame(rows)
    cath_df.to_csv(f"{args.out_prefix}_per_protein.csv",index=False)

    if args.input_table:
        tbl=pd.read_csv(args.input_table)
        tbl["_pdb_id"]=tbl[args.protein_col].astype(str).str.extract(r"^([0-9A-Za-z]{4})",expand=False).str.lower()
        annotated=tbl.merge(cath_df,on="pdb_id",how="left")
        annotated.to_csv(f"{args.out_prefix}_annotated.csv",index=False)

if __name__=="__main__":
    main()
