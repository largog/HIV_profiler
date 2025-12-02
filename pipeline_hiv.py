import os
import pandas as pd
from pathlib import Path
import subprocess
import json
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem import Crippen
from rdkit.Chem import Draw


# ============================================================
#                MÓDULO 1 — BLAST + IDENTIFICACIÓN
# ============================================================

def run_blast(query_fasta, db_path, out_tsv="blast_output.tsv"):
    cmd = [
        "blastn",
        "-query", str(query_fasta),
        "-db", str(db_path),
        "-outfmt", (
            "6 qseqid sseqid pident length qlen slen "
            "qstart qend sstart send evalue bitscore"
        ),
        "-max_target_seqs", "5",
        "-out", out_tsv
    ]
    subprocess.run(cmd, check=True)
    return Path(out_tsv)


def extract_subtype(header: str):
    header = header.lstrip(">")
    parts = header.split(".")
    if len(parts) >= 2:
        return parts[1]   # ejemplo: 01_AE
    return "Unknown"


def parse_top_hit(blast_tsv):
    cols = [
        "qseqid", "sseqid", "pident", "length", "qlen", "slen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore"
    ]

    df = pd.read_csv(blast_tsv, sep="\t", names=cols)
    if df.empty:
        return None

    top = df.sort_values("bitscore", ascending=False).iloc[0]
    coverage = top["length"] / top["slen"]
    subtype = extract_subtype(top["sseqid"])

    return {
        "query_id": top["qseqid"],
        "reference_id": top["sseqid"],
        "reference_length":int(top["slen"]),
        "subtype": subtype,
        "percent_identity": float(round(top["pident"], 3)),
        "alignment_length": int(top["length"]),
        "query_length": int(top["qlen"]),
        "coverage": float(round(coverage, 3)),
        "evalue": float(top["evalue"]),
        "bitscore": float(top["bitscore"]),
        "S_start": int(top["sstart"]),
        "S_end": int(top["send"]),
    }


def identify_variant(query_fasta, db_path):
    blast_output = run_blast(query_fasta, db_path)
    return parse_top_hit(blast_output)



# ============================================================
#            MÓDULO 2 — DATOS FARMACOLÓGICOS
# ============================================================

def load_subtype_pharma_data(blast_result, data_folder="."):
    subtype = blast_result["subtype"]
    filename = f"{subtype}.tsv"
    filepath = os.path.join(data_folder, filename)

    if not os.path.exists(filepath):
        raise FileNotFoundError(f"No existe un archivo para el subtipo: {filepath}")

    df = pd.read_csv(filepath, sep="\t")

    important_cols = [
        "BindingDB Reactant_set_id",
        "Ligand SMILES",
        "BindingDB Ligand Name",
        "Target Name",
        "Ki (nM)",
        "IC50 (nM)",
        "Kd (nM)",
        "EC50 (nM)",
    ]

    important_cols = [c for c in important_cols if c in df.columns]

    return df[important_cols]



# ============================================================
#            MÓDULO 3 — TOP 3 MEJORES FÁRMACOS
# ============================================================

def get_top_drugs(pharma_df, top_n=3):
    """
    Selecciona los mejores compuestos según:
      1. pChEMBL (mayor = mejor)
      2. Ki (menor = mejor)
      3. IC50 (menor = mejor)
    """

    df = pharma_df.copy()

    numeric_cols = ["pChEMBL Value", "Ki (nM)", "IC50 (nM)"]
    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    sort_order = []

    if "pChEMBL Value" in df.columns:
        sort_order.append(("pChEMBL Value", False))  # mayor mejor

    if "Ki (nM)" in df.columns:
        sort_order.append(("Ki (nM)", True))  # menor mejor

    if "IC50 (nM)" in df.columns:
        sort_order.append(("IC50 (nM)", True))  # menor mejor

    if not sort_order:
        raise ValueError("No existen columnas pChEMBL, Ki o IC50 para ordenar.")

    df_sorted = df.sort_values(
        by=[c for c, _ in sort_order],
        ascending=[a for _, a in sort_order]
    )

    return df_sorted.head(top_n)



# ============================================================
#            MÓDULO 4 — ANÁLISIS RDKit (TOP 3)
# ============================================================

def analyze_top3_with_rdkit(top3, output_folder="results"):
    results = []

    os.makedirs(output_folder, exist_ok=True)

    for i, row in top3.iterrows():
        smiles = row.get("Ligand SMILES")
        raw = row.get("BindingDB Ligand Name")

        commercial_name = "Unknown"
        iupac_name = None

        if isinstance(raw, str):
            parts = [p.strip() for p in raw.split("::")]

            # --- Detectar nombre comercial (más confiable) ---
            for p in reversed(parts):
                if p.startswith("CHEMBL"):
                    continue
                if any(char.isdigit() for char in p):  # evitar códigos ABT-538, CGP 73547, etc.
                    continue
                if p.isupper() and len(p) <= 4:  # LPV, RTV, etc
                    continue

                if p[0].isupper() and p[1:].islower():  # patrón de nombre de fármaco
                    commercial_name = p
                    break

            # --- Detectar IUPAC real ---
            # regla: el PRIMER fragmento químico largo es el IUPAC
            for p in parts:
                if (
                    len(p) > 25 and
                    any(x in p for x in ["(", ")", "-", "yl", "oxy", "methyl", "phenyl", "hydroxy"])
                ):
                    iupac_name = p
                    break

        if iupac_name is None:
            iupac_name = f"Mol_{i+1}"

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue

        image_path = os.path.join(output_folder, f"top3_mol_{i+1}.png")
        Draw.MolToFile(mol, image_path, size=(400, 400))

        props = {
            "Index": i + 1,
            "Commercial Name": commercial_name,
            "IUPAC Name": iupac_name,
            "SMILES": smiles,
            "MW": Descriptors.MolWt(mol),
            "LogP": Crippen.MolLogP(mol),
            "HBA": rdMolDescriptors.CalcNumHBA(mol),
            "HBD": rdMolDescriptors.CalcNumHBD(mol),
            "TPSA": rdMolDescriptors.CalcTPSA(mol),
            "Rotatable Bonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
            "Rings": rdMolDescriptors.CalcNumRings(mol),
            "Formula": rdMolDescriptors.CalcMolFormula(mol),
            "pChEMBL": row.get("pChEMBL Value"),
            "Ki (nM)": row.get("Ki (nM)"),
            "IC50 (nM)": row.get("IC50 (nM)"),
            "Image": image_path
        }

        results.append(props)

    return pd.DataFrame(results)


# ============================================================
#            MÓDULO 5 — GUARDAR RESULTADOS
# ============================================================

def save_results(blast_info, pharma_df, top3, output_folder="results"):
    os.makedirs(output_folder, exist_ok=True)

    with open(os.path.join(output_folder, "blast_result.json"), "w") as f:
        json.dump(blast_info, f, indent=4)

    pharma_df.to_csv(os.path.join(output_folder, "pharma_data.csv"), index=False)
    top3.to_csv(os.path.join(output_folder, "top3_drugs.csv"), index=False)

    print(f"\n📁 Resultados guardados en: {output_folder}")



# ============================================================
#                    PIPELINE COMPLETO
# ============================================================

def full_pipeline(query_fasta, db_path, data_folder="."):
    print("\nEjecutando BLAST...")
    blast_info = identify_variant(query_fasta, db_path)
    print("✔ Subtipo identificado:", blast_info["subtype"])

    print("\nCargando archivo farmacológico...")
    pharma_df = load_subtype_pharma_data(blast_info, data_folder)

    print("\nSeleccionando Top 3...")
    top3 = get_top_drugs(pharma_df)

    print("\nRDKit...")
    rdkit_info = analyze_top3_with_rdkit(top3)

    print("\nGuardando imágenes...")
    save_results(blast_info, pharma_df, top3)

    return blast_info, pharma_df, top3, rdkit_info



# ============================================================
#                    EJECUCIÓN LOCAL
# ============================================================

if __name__ == "__main__":
    QUERY = "data_test/B1.fasta"
    DATABASE = "database/HIV1_REF"

    blast_info, pharma_df, top3, rdkit_info = full_pipeline(
        QUERY, DATABASE, data_folder="."
    )

    print("\n==== BLAST ====\n", blast_info)
    print("\n==== PHARMA ====\n", pharma_df.head())
    print("\n==== TOP3 ====\n", top3)
    print("\n==== RDKIT ====\n", rdkit_info)




