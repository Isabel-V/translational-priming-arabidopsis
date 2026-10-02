"""Recompute every number quoted in the manuscript from the raw tables.

Writes results/manuscript_numbers.tsv (claim, value) and results/master_gene_table.tsv.
"""
import numpy as np
import pandas as pd
from te_common import master_table, RESULTS, TE_HIGH

RESULTS.mkdir(exist_ok=True)
m = master_table()
m.to_csv(RESULTS / "master_gene_table.tsv", sep="\t", index=False)
N = []
def rec(k, v): N.append((k, v)); print(f"{k:<75} {v}")

rec("gene models", len(m))
rna_det = m.rna_detected_root | m.rna_detected_shoot
ribo_det = m.ribo_detected_root | m.ribo_detected_shoot
rec("RNA-seq detected (RPKM>0, root or shoot)", int(rna_det.sum()))
rec("Ribo-seq detected", int(ribo_det.sum()))
rec("detected in both RNA-seq and Ribo-seq", int((rna_det & ribo_det).sum()))
rec("RNA detected root / shoot", (int(m.rna_detected_root.sum()), int(m.rna_detected_shoot.sum())))

# TE range (non-zero TE only)
for t in ("root", "shoot"):
    te = m[f"te_{t}"].where(m[f"te_{t}"] > 0)
    i, j = te.idxmin(), te.idxmax()
    rec(f"TE range {t}", f"{te.min():.4g} ({m.gene_model[i]}) - {te.max():.6g} ({m.gene_model[j]}); fold {te.max()/te.min():.3g}")
allte = pd.concat([m.te_root, m.te_shoot]); allte = allte[allte > 0]
rec("TE range both tissues combined, fold", f"{allte.max()/allte.min():.3g}")
for g in ("AT4G21400", "AT3G55860"):
    r = m[m.agi == g].iloc[0]
    rec(f"TE {g} root/shoot", (round(r.te_root, 4), None if pd.isna(r.te_shoot) else round(r.te_shoot, 4)))

# TE classes (union)
rec("high-TE any tissue", int(m.high_te_any.sum()))
rec("mid-TE any tissue", int(m.mid_te_any.sum()))
rec("low-TE any tissue", int(m.low_te_any.sum()))
lo_real = ((m.te_class_root == "low-TE") & ~m.low_te_zero_ribo_root) | ((m.te_class_shoot == "low-TE") & ~m.low_te_zero_ribo_shoot)
rec("low-TE any tissue, excluding TE=0 (no Ribo-seq reads)", int(lo_real.sum()))
rec("low-TE only because Ribo-seq = 0 in that tissue", int(m.low_te_any.sum() - lo_real.sum()))
rec("high-TE both tissues", int(m.high_te_both.sum()))
rec("baseline high-TE fraction (/27,655)", round(m.high_te_any.mean() * 100, 2))

# RNA-low population
low_rna = ((m.rna_rpkm_root > 0) & (m.rna_rpkm_root < 1)) | ((m.rna_rpkm_shoot > 0) & (m.rna_rpkm_shoot < 1))
rec("genes with 0<RPKM<1 in >=1 tissue", int(low_rna.sum()))
for t in ("root", "shoot"):
    s = m[(m[f"rna_rpkm_{t}"] > 0) & (m[f"rna_rpkm_{t}"] < 1)][f"rna_reads_{t}"]
    rec(f"RNA reads for 0<RPKM<1 genes, {t}: median [IQR]", f"{s.median():.0f} [{s.quantile(.25):.0f}-{s.quantile(.75):.0f}]")

# primed
rec("primed root / shoot", (int(m.primed_root.sum()), int(m.primed_shoot.sum())))
rec("primed union", int(m.primed_any.sum()))
rec("primed both tissues", int(m.primed_both.sum()))
rec("primed both, nuclear / chloroplast / mito", m[m.primed_both].organelle.value_counts().to_dict())
rec("primed union nuclear", int((m.primed_any & (m.organelle == "nuclear")).sum()))
rec("primed union as % of high-TE", round(100 * m.primed_any.sum() / m.high_te_any.sum(), 1))
canon_only = m.high_te_any & ~m.primed_any
rec("high-TE never primed (canonical in all high-TE tissues)", int(canon_only.sum()))
# tissue-split: primed in one tissue, high-TE with RNA>=0 in the other
split = m.primed_any & (((m.te_class_root == "high-TE") & ~m.primed_root) | ((m.te_class_shoot == "high-TE") & ~m.primed_shoot))
rec("primed genes that are canonical high-TE in the other tissue", int(split.sum()))
rec("primed-only (no canonical high-TE tissue)", int((m.primed_any & ~split).sum()))
rec("mean log2TE primed root / shoot", (round(m.loc[m.primed_root, 'log2_te_root'].mean(), 2), round(m.loc[m.primed_shoot, 'log2_te_shoot'].mean(), 2)))

# proteomics
rec("protein loci detected (High+Medium) root/shoot/any", (int(m.protein_root.sum()), int(m.protein_shoot.sum()), int(m.protein_any.sum())))
p = m[m.primed_any]
rec("primed with protein evidence", f"{int(p.protein_any.sum())}/{len(p)} = {100*p.protein_any.mean():.1f}%")
pn = p[p.organelle == "nuclear"]
rec("nuclear primed with protein evidence", f"{int(pn.protein_any.sum())}/{len(pn)} = {100*pn.protein_any.mean():.1f}%")
pb = m[m.primed_both]
rec("primed-both with protein evidence", f"{int(pb.protein_any.sum())}/{len(pb)} = {100*pb.protein_any.mean():.1f}%")
for lab, mask in (("all detected genes (RNA+Ribo)", rna_det & ribo_det), ("canonical high-TE (never primed)", canon_only),
                  ("mid-TE any", m.mid_te_any), ("low-TE any (with Ribo reads)", lo_real)):
    s = m[mask]
    rec(f"protein evidence, {lab}", f"{int(s.protein_any.sum())}/{len(s)} = {100*s.protein_any.mean():.1f}%")

# canonical reference sets for sequence analysis
rec("seq sets (log2RNA>=2 both; TE both)", m.seq_set.value_counts().to_dict())
expr = (m.log2_rna_root >= 2) & (m.log2_rna_shoot >= 2)
rec("genes with log2RNA>=2 in both tissues", int(expr.sum()))

# micropeptides among 41 nuclear primed-both
s41 = m[m.primed_both & (m.organelle == "nuclear")].copy()
s41["mean_log2_te"] = s41[["log2_te_root", "log2_te_shoot"]].mean(axis=1)
s41 = s41.sort_values("mean_log2_te", ascending=False)
rec("41 nuclear shared: top 6 by mean TE", s41[["gene_model", "predicted_protein_length_aa", "mean_log2_te"]].head(6).round(2).values.tolist())
rec("41 nuclear shared: n with predicted protein <200 aa", int((s41.predicted_protein_length_aa < 200).sum()))

# iron genes
iron = {"IRT1": "AT4G19690", "IRT2": "AT4G19680", "FRO2": "AT1G01580", "VTL1": "AT1G21140", "VTL2": "AT1G76800",
        "AT3G43630 (VIT family)": "AT3G43630", "MTP8": "AT3G58060", "S8H": "AT3G12900", "F6'H1": "AT3G13610",
        "CYP82C4": "AT4G31940", "CYP82C3": "AT4G31950", "CYP82C2": "AT4G31970", "NAS2": "AT5G56080", "NAS3": "AT1G09240",
        "FIT": "AT2G28160", "bHLH38": "AT3G56970", "bHLH39": "AT3G56980", "bHLH100": "AT2G41240", "bHLH101": "AT5G04150",
        "PYE": "AT3G47640", "BTS": "AT3G18290", "PDR9": "AT3G53480", "MYB72": "AT1G56160"}
for k, g in iron.items():
    r = m[m.agi == g].iloc[0]
    rec(f"iron {k} {g}: RNA r/s, TE r/s, primed r/s, protein",
        (round(r.log2_rna_root, 2), round(r.log2_rna_shoot, 2), round(r.log2_te_root, 2), round(r.log2_te_shoot, 2),
         bool(r.primed_root), bool(r.primed_shoot), bool(r.protein_any)))

# S8H rank
for t in ("root", "shoot"):
    x = m[f"log2_te_{t}"].dropna()
    v = m.loc[m.agi == "AT3G12900", f"log2_te_{t}"].iloc[0]
    rec(f"S8H log2TE rank {t}", f"{int((x > v).sum()) + 1} of {len(x)} (top {100*(x > v).mean():.2f}%)")

pd.DataFrame(N, columns=["quantity", "value"]).to_csv(RESULTS / "manuscript_numbers.tsv", sep="\t", index=False)
