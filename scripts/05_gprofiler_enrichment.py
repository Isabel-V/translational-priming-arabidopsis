"""GO and KEGG enrichment (g:Profiler g:GOSt) for the gene sets cited in the manuscript.

Step 1 (always, offline): writes the gene lists from results/master_gene_table.tsv (script 01) to
    data/gprofiler/gene_lists/ and checks them against the lists used for the published run.
Step 2 (only with --run; needs internet access to biit.cs.ut.ee):
    pip install gprofiler-official
    python 05_gprofiler_enrichment.py --run
    -> data/gprofiler/results/gprofiler_<set>_{all,significant}.csv, gprofiler_summary.tsv, run_info.txt

The results distributed with the repository were produced with g:Profiler version e114_eg62_p19_27110d83
(run 2026-09-28). Parameters (identical for all sets): organism athaliana; GO:BP, GO:MF, GO:CC, KEGG;
Benjamini-Hochberg FDR < 0.05 (g:SCS not used); electronic GO annotations (IEA) excluded;
background = annotated genome. The primed set is additionally tested against all nuclear genes detected by
RNA-seq (custom background, sensitivity analysis).
"""
import sys
import datetime
import pandas as pd
from te_common import ROOT, RESULTS

GP = ROOT / "data" / "gprofiler"
LISTS, OUT = GP / "gene_lists", GP / "results"
LISTS.mkdir(parents=True, exist_ok=True); OUT.mkdir(exist_ok=True)

m = pd.read_csv(RESULTS / "master_gene_table.tsv", sep="\t")
nuc = m.organelle == "nuclear"
canon = m.high_te_any & ~m.primed_any
SETS = {
    "primed_nuclear_701": m.primed_any & nuc,             # translationally primed in >= 1 tissue
    "primed_both_nuclear_41": m.primed_both & nuc,        # primed in roots and shoots
    "highTE_both_nuclear_93": m.high_te_both & nuc,       # high-TE in both tissues
    "highTE_any_nuclear_855": m.high_te_any & nuc,        # all nuclear high-TE genes
    "canonical_highTE_nuclear_154": canon & nuc,          # high-TE, never primed
    "background_detected_nuclear_24755": (m.rna_detected_root | m.rna_detected_shoot) & nuc,
}
lists = {}
for name, mask in SETS.items():
    genes = sorted(m[mask].agi.unique())
    f = LISTS / f"{name}.txt"
    if f.exists():
        old = sorted(g for g in f.read_text().split() if g)
        status = "identical to distributed list" if old == genes else "DIFFERS from distributed list"
    else:
        status = "new"
    f.write_text("\n".join(genes) + "\n")
    lists[name] = genes
    print(f"{name:36s} {len(genes):6d} genes  ({status})")

if "--run" not in sys.argv:
    print("\nGene lists written. Add --run to query g:Profiler (internet required).")
    sys.exit(0)

from gprofiler import GProfiler
PARAMS = dict(organism="athaliana", sources=["GO:BP", "GO:MF", "GO:CC", "KEGG"], user_threshold=0.05,
              significance_threshold_method="fdr", no_iea=True, no_evidences=False)
gp = GProfiler(return_dataframe=True, user_agent="te_priming_manuscript")
background = lists["background_detected_nuclear_24755"]
runs = [(s, lists[s], None) for s in list(SETS)[:-1]] + \
       [("primed_nuclear_701_bgDetected", lists["primed_nuclear_701"], background)]
summary, version = [], None
for name, genes, bg in runs:
    kw = dict(PARAMS)
    if bg is not None:
        kw["background"] = bg
    res = gp.profile(query=genes, **kw)
    version = version or (gp.meta or {}).get("version")
    res = res.sort_values("p_value") if len(res) else res
    res.to_csv(OUT / f"gprofiler_{name}_all.csv", index=False)
    sig = res[res["significant"]] if len(res) else res
    sig.to_csv(OUT / f"gprofiler_{name}_significant.csv", index=False)
    qs = res.groupby("source")["query_size"].first().to_dict() if len(res) else {}
    summary.append(dict(set=name, n_genes=len(genes), n_significant=len(sig),
                        significant_by_source=sig["source"].value_counts().to_dict() if len(sig) else {},
                        annotated_genes_per_source=qs, top_terms="; ".join(sig["name"].head(10)) if len(sig) else ""))
    print(f"{name:34s} {len(genes):5d} genes  {len(sig):3d} significant terms")
pd.DataFrame(summary).to_csv(OUT / "gprofiler_summary.tsv", sep="\t", index=False)
info = (f"date: {datetime.date.today()}\ng:Profiler version: {version}\nparameters: {PARAMS}\n"
        f"custom background (sensitivity run): {len(background)} nuclear genes detected by RNA-seq\n")
(OUT / "run_info.txt").write_text(info)
print("\n" + info)
