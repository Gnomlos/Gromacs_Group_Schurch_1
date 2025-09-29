import os
import time
import requests
import pandas as pd

# This script allows to extract the pdb files, pdbID and UniprotID out of the tsv file form the HPA librairy. 


# ==== CONFIG ====
HPA_FILE = "/home/gadolinium/Documents/Protein pipeline/sa_location_Secreted.tsv"  # your HPA file -- change if you want to change studied sector.
UNIPROT_COLUMN = "Uniprot"  # change to the actual column name in your TSV
OUTPUT_DIR = "/home/gadolinium/Documents/Protein pipeline/pdb_structures"  # output folder to save downloaded structures
FILE_FORMAT = "pdb"  # "pdb" or "cif"
MAPPING_FILE = os.path.join(OUTPUT_DIR, "uniprot_pdb_cath_mapping.csv")

RCSB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
UNIPROT_ENTRY_URL = "https://rest.uniprot.org/uniprotkb/{acc}.json"

# Make output folder
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Load HPA TSV
df = pd.read_csv(HPA_FILE, sep="\t")

# Extract UniProt IDs
uniprot_ids = (
    df[UNIPROT_COLUMN]
    .dropna()
    .astype(str)
    .str.strip()
    .str.split(r"[;,\s]+")   # handle multi-IDs in one cell if present
    .explode()
    .dropna()
    .unique()
)

def fetch_pdb_ids(uniprot_id):
    """Query RCSB PDB for all structures of a given UniProt ID."""
    query = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                "operator": "exact_match",
                "value": uniprot_id
            }
        },
        "return_type": "entry"
    }
    try:
        r = requests.post(RCSB_SEARCH_URL, json=query, timeout=30)
        if r.status_code != 200:
            return []
        results = r.json().get("result_set", [])
        return [res["identifier"] for res in results]
    except requests.RequestException:
        return []

def download_structure(pdb_id, fmt="pdb"):
    """Download a structure (PDB or CIF) by its ID."""
    url = f"https://files.rcsb.org/download/{pdb_id}.{fmt}"
    try:
        r = requests.get(url, timeout=60)
        if r.status_code == 200:
            file_path = os.path.join(OUTPUT_DIR, f"{pdb_id}.{fmt}")
            # text write is fine for .pdb/.cif
            with open(file_path, "w", encoding="utf-8") as f:
                f.write(r.text)
            print(f"Downloaded {pdb_id}.{fmt}")
        else:
            print(f"Failed to download {pdb_id}.{fmt} (HTTP {r.status_code})")
    except requests.RequestException as e:
        print(f"Failed to download {pdb_id}.{fmt}: {e}")

def uniprot_to_cath(acc, retries=3, pause=0.5):
    """
    Return a set of CATH superfamily IDs for a UniProt accession,
    based on UniProt cross-references (Gene3D / CATH id).
    """
    url = UNIPROT_ENTRY_URL.format(acc=acc)
    for _ in range(retries):
        try:
            r = requests.get(url, headers={"User-Agent":"uniprot-cath-mapper"}, timeout=30)
            if r.status_code == 200:
                data = r.json()
                cath_ids = set()
                for xref in data.get("uniProtKBCrossReferences", []):
                    if xref.get("database") == "Gene3D":
                        props = {p.get("key"): p.get("value") for p in xref.get("properties", [])}
                        # Prefer explicit "CATH id" if present; else fall back to xref id
                        cid = props.get("CATH id") or xref.get("id")
                        if cid:
                            cath_ids.add(cid)
                return cath_ids
        except requests.RequestException:
            pass
        time.sleep(pause)
    return set()

def save_mapping(mapping, filepath):
    """Save UniProt–PDB–CATH mapping to a CSV file."""
    df_map = pd.DataFrame(mapping, columns=["UniProt_ID", "PDB_ID", "CATH_Superfamily_IDs"])
    df_map.to_csv(filepath, index=False)
    print(f"Mapping saved to {filepath}")

# --- Main ---
mapping_records = []

for i, uid in enumerate(uniprot_ids, 1):
    uid = uid.strip()
    if not uid:
        continue

    pdb_ids = fetch_pdb_ids(uid)

    # Fetch CATH IDs once per UniProt accession (reuse for all its PDBs)
    cath_ids = uniprot_to_cath(uid)
    cath_joined = ";".join(sorted(cath_ids)) if cath_ids else ""

    for pdb_id in pdb_ids:
        download_structure(pdb_id, FILE_FORMAT)
        mapping_records.append((uid, pdb_id, cath_joined))

    if i % 50 == 0:
        print(f"...processed {i}/{len(uniprot_ids)} UniProt accessions")

# Save mapping at the end
if mapping_records:
    save_mapping(mapping_records, MAPPING_FILE)
else:
    print("No mappings found. Nothing to save.")
