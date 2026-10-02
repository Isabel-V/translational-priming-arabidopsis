"""Shared definitions for the translational-priming analysis.

All thresholds used in the manuscript are defined here once, so that every
script (tables, figures, verification) applies exactly the same criteria.
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
SI = ROOT / "supporting_information"

TE_FILE = RAW / "TranslationEfficiency2.xlsx"          # Wen-Dar Lin, processed RNA-seq/Ribo-seq
PROT_FILE = RAW / "Jain_Schmidt_2023_mmc3.xlsx"         # Jain & Schmidt (2023) Mol Cell Proteomics, Suppl. dataset
ANNOT_FILE = RAW / "TE_full_annotation.xlsx"            # Araport11 descriptions + functional categories
MOTIF_FILE = RAW / "motif_analysis_Supplementary_Data_5.xlsx"

# ---- thresholds (log2 scale) -------------------------------------------
TE_HIGH = 1.5            # high-TE:  log2TE >  1.5
TE_LOW = -1.5            # low-TE:   log2TE < -1.5
RNA_PRIMED = 0.0         # translationally primed: log2 RNA-seq RPKM < 0 (RPKM < 1)
RNA_CANONICAL_MOTIF = 2.0  # canonical reference sets for sequence analysis: log2 RPKM >= 2 in both tissues
ZERO_CODE = -99          # Wen-Dar Lin's files code log2(0) as -100

TISSUES = ("root", "shoot")


def _log2(x):
    x = np.asarray(x, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(x > 0, np.log2(x), np.nan)


def load_te_table() -> pd.DataFrame:
    """Per-gene RPKM, read counts, log2 values and TE for roots and shoots."""
    out = None
    for tissue, sheet in (("root", "root-ctr-0h (2)"), ("shoot", "shoot-ctr-0h (2)")):
        d = pd.read_excel(TE_FILE, sheet_name=sheet, usecols=range(17))
        d = d.rename(columns={
            "#geneID": "gene_model", "#read": f"rna_reads_{tissue}",
            "#ribo": f"ribo_reads_5utr_{tissue}", "#ribo.1": f"ribo_reads_cds_{tissue}",
            "#ribo.2": f"ribo_reads_3utr_{tissue}", "rnaseqRPKM": f"rna_rpkm_{tissue}",
            "riboRPKM": f"ribo_rpkm_{tissue}", "cdsLen": "cds_length_nt",
        })
        keep = ["gene_model", "cds_length_nt", f"rna_reads_{tissue}", f"ribo_reads_5utr_{tissue}",
                f"ribo_reads_cds_{tissue}", f"ribo_reads_3utr_{tissue}", f"rna_rpkm_{tissue}", f"ribo_rpkm_{tissue}"]
        d = d[keep]
        if out is None:
            out = d
        else:
            out = out.merge(d.drop(columns="cds_length_nt"), on="gene_model", how="outer")
    out["agi"] = out["gene_model"].str.split(".").str[0]
    for t in TISSUES:
        out[f"ribo_reads_{t}"] = out[[f"ribo_reads_5utr_{t}", f"ribo_reads_cds_{t}", f"ribo_reads_3utr_{t}"]].sum(axis=1)
        out[f"log2_rna_{t}"] = _log2(out[f"rna_rpkm_{t}"])
        out[f"log2_ribo_{t}"] = _log2(out[f"ribo_rpkm_{t}"])
        rna, ribo = out[f"rna_rpkm_{t}"], out[f"ribo_rpkm_{t}"]
        te = np.where(rna > 0, ribo / rna.where(rna > 0), np.nan)
        out[f"te_{t}"] = te
        out[f"log2_te_{t}"] = _log2(te)                      # NaN when ribo = 0 or rna = 0
        out[f"rna_detected_{t}"] = rna > 0
        out[f"ribo_detected_{t}"] = ribo > 0
        # TE class as in Wen-Dar Lin's classification (TE = 0 counted as low-TE)
        cls = pd.Series(pd.NA, index=out.index, dtype="object")
        l2 = out[f"log2_te_{t}"]
        cls[(rna > 0) & (ribo == 0)] = "low-TE"
        cls[l2 < TE_LOW] = "low-TE"
        cls[(l2 >= TE_LOW) & (l2 <= TE_HIGH)] = "mid-TE"
        cls[l2 > TE_HIGH] = "high-TE"
        out[f"te_class_{t}"] = cls
        out[f"low_te_zero_ribo_{t}"] = (rna > 0) & (ribo == 0)
        out[f"primed_{t}"] = (out[f"log2_rna_{t}"] < RNA_PRIMED) & (l2 > TE_HIGH)
    out["organelle"] = np.select([out.agi.str.startswith("ATCG"), out.agi.str.startswith("ATMG")],
                                 ["chloroplast", "mitochondrion"], "nuclear")
    out["predicted_protein_length_aa"] = (out["cds_length_nt"] // 3 - 1).astype("Int64")
    return out


def load_protein_evidence() -> pd.DataFrame:
    """AGI loci detected (High or Medium FDR confidence) in the Jain & Schmidt (2023) proteomes."""
    rows = []
    for tissue, sheet in (("root", "Root Proteome"), ("shoot", "Shoot Proteome")):
        d = pd.read_excel(PROT_FILE, sheet_name=sheet, header=1)
        d = d[d["Protein FDR Confidence: Combined"].isin(["High", "Medium"])]
        loci = d["Accession"].astype(str).str.upper().str.split(".").str[0]
        loci = loci[loci.str.match(r"^AT[1-5CM]G\d{5}$")]
        rows.append(pd.DataFrame({"agi": loci.unique(), f"protein_{tissue}": True}))
    p = rows[0].merge(rows[1], on="agi", how="outer")
    p = p.fillna(False).infer_objects(copy=False)
    p["protein_any"] = p.protein_root | p.protein_shoot
    return p


def load_annotation() -> pd.DataFrame:
    frames = [pd.read_excel(ANNOT_FILE, sheet_name=s, header=3) for s in ("HighTE_All", "MidTE_All", "LowTE_All")]
    a = pd.concat(frames).drop_duplicates("Gene ID")
    return a.rename(columns={"Gene ID": "gene_model", "Description": "description",
                             "Functional category": "functional_category"})[["gene_model", "description", "functional_category"]]


def master_table() -> pd.DataFrame:
    te = load_te_table()
    prot = load_protein_evidence()
    ann = load_annotation()
    m = te.merge(prot, on="agi", how="left").merge(ann, on="gene_model", how="left")
    for c in ("protein_root", "protein_shoot", "protein_any"):
        m[c] = m[c].fillna(False).astype(bool)
    m["primed_any"] = m.primed_root | m.primed_shoot
    m["primed_both"] = m.primed_root & m.primed_shoot
    m["high_te_any"] = (m.te_class_root == "high-TE") | (m.te_class_shoot == "high-TE")
    m["mid_te_any"] = (m.te_class_root == "mid-TE") | (m.te_class_shoot == "mid-TE")
    m["low_te_any"] = (m.te_class_root == "low-TE") | (m.te_class_shoot == "low-TE")
    m["high_te_both"] = (m.te_class_root == "high-TE") & (m.te_class_shoot == "high-TE")
    # canonical sequence-analysis reference sets (both tissues, log2 RNA >= 2)
    expr = (m.log2_rna_root >= RNA_CANONICAL_MOTIF) & (m.log2_rna_shoot >= RNA_CANONICAL_MOTIF)
    lr, ls = m.log2_te_root, m.log2_te_shoot
    lr0 = lr.where(~m.low_te_zero_ribo_root, -np.inf)
    ls0 = ls.where(~m.low_te_zero_ribo_shoot, -np.inf)
    hi = (lr0 > TE_HIGH) & (ls0 > TE_HIGH)
    lo = (lr0 < TE_LOW) & (ls0 < TE_LOW)
    mid = (lr0 <= TE_HIGH) & (ls0 <= TE_HIGH) & ~lo
    m["seq_set"] = np.select([expr & hi, expr & lo, expr & mid], ["high-TE", "low-TE", "low/mid-TE"], "")
    m.loc[~expr, "seq_set"] = ""
    return m
