#!/usr/bin/env python3
# This script is used to to insert the Gadolinium ion in the pocket, in the pdb file so that the simulation can take place.

import csv, sys, re, pathlib

pdb_in   = pathlib.Path(sys.argv[1])     # protein.pdb (pre-protonated is best)
csv_in   = pathlib.Path(sys.argv[2])     # pockets.csv with columns x,y,z or center_x,center_y,center_z
tag      = sys.argv[3] if len(sys.argv) > 3 else "POCKET"

pdb = pdb_in.read_text().splitlines()
# find last atom serial & last res seq for robust numbering
serial = max([int(l[6:11]) for l in pdb if l.startswith(("ATOM  ","HETATM"))] or [0])
resseq = max([int(re.sub(r"\D","",l[22:26]) or 0) for l in pdb if l.startswith(("ATOM  ","HETATM"))] or [0])

def hetatm_line(sn, resn, x,y,z):
    # PDB: columns must be aligned; residue name 'GD3', atom name 'GD', element 'Gd'
    return f"HETATM{sn:5d} {'GD':>4} {'GD3':>3} {'X':1}{resn:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          Gd"

with csv_in.open() as f:
    rdr = csv.DictReader(f)
    cols = rdr.fieldnames
    # flexible column names
    if {"x","y","z"}.issubset(set(c.lower() for c in cols)):
        get = lambda r: (float(r["x"]), float(r["y"]), float(r["z"]))
    else:
        get = lambda r: (float(r.get("center_x") or r.get("cx") or r.get("X")),
                         float(r.get("center_y") or r.get("cy") or r.get("Y")),
                         float(r.get("center_z") or r.get("cz") or r.get("Z")))
    for i,row in enumerate(rdr, start=1):
        x,y,z = get(row)
        serial += 1; resseq += 1
        out = []
        out.extend(pdb)
        out.append(hetatm_line(serial, resseq, x,y,z))
        out.append("TER")
        out.append("END")
        out_pdb = pdb_in.with_name(f"{pdb_in.stem}_{tag}{i:03d}_GD3.pdb")
        out_pdb.write_text("\n".join(out))
        print(f"Wrote {out_pdb}")


# bash-use: python insert_gd.py protein_pH7.pdb pockets.csv
