"""Re-analysis of the CTATAWAWA-like motif (MEME motif 6) in start-codon-proximal windows.

Inputs (data/motif/):
  high-TE_200bp.fasta, low-TE_200bp.fasta, low-mid-TE_200bp.fasta  (-150..+50 relative to the A of the AUG;
                                                                      A of AUG at 0-based index 150)
  CTATAWAWA_motif.meme   (MEME letter-probability matrix, discovered in the 27 high-TE sequences)

What it does
  1. Scores every window with the motif PWM (log-odds vs. the nucleotide composition of all windows,
     FIMO-style pseudocount), sense strand; a site = window position with PWM P < 1e-4.
  2. Tests enrichment (two-sided Fisher's exact test) in
       - the discovery set (high-TE, n = 27; circular, the motif was built from these sequences),
       - two independent validation sets not used for discovery: 1.0 < log2TE <= 1.5 and 0.5 < log2TE <= 1.5
         in both tissues (sequences taken from the low/mid-TE file),
     against low-TE genes and against the remaining low/mid-TE genes.
  3. Compares each set with first-order Markov (dinucleotide) shuffles of its own sequences (SEA-like control).
  4. Locates every site relative to the AUG and to the annotated transcription start site (5'UTR length from
     Araport11 models in TranslationEfficiency2.xlsx) -> within the 5'UTR (transcribed) or upstream (promoter).
  5. Scans the -150..-1 window of the 701 nuclear translationally primed genes (data/motif/primed_701_upstream150.fasta)
     and compares it with the same window of the sets above; sites are located relative to the annotated TSS.
Outputs: results/motif_enrichment.tsv, motif_sites.tsv, motif_primed_upstream150.tsv,
         figures/FigureS1_motif_reanalysis.{png,pdf}
"""
import re
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from te_common import ROOT, RESULTS, FIGURES, TE_FILE, master_table

MOTIF = ROOT / "data" / "motif"
RNG = np.random.default_rng(1)
P_SITE = 1e-4
ATG = 150
BASES = "ACGT"


def read_fasta(f):
    d, k = {}, None
    for line in open(f):
        line = line.strip()
        if line.startswith(">"):
            k = line[1:]; d[k] = []
        elif k:
            d[k].append(line.upper())
    return {k: "".join(v) for k, v in d.items()}


def read_meme(f):
    rows, on = [], False
    for line in open(f):
        if line.startswith("letter-probability"):
            on = True; nsites = int(re.search(r"nsites= *(\d+)", line).group(1)); continue
        if on:
            v = line.split()
            if len(v) == 4:
                rows.append([float(x) for x in v])
            elif rows:
                break
    return np.array(rows), nsites


seqs = {"high": read_fasta(MOTIF / "high-TE_200bp.fasta"), "low": read_fasta(MOTIF / "low-TE_200bp.fasta"),
        "lowmid": read_fasta(MOTIF / "low-mid-TE_200bp.fasta")}
allseq = "".join(s for d in seqs.values() for s in d.values())
bg = np.array([allseq.count(b) for b in BASES], float); bg /= bg.sum()
P, nsites = read_meme(MOTIF / "CTATAWAWA_motif.meme")
W = len(P)
Pp = (P * nsites + 0.1 * bg) / (nsites + 0.1)            # FIMO-style pseudocount
LO = np.log2(Pp / bg)                                       # W x 4 log-odds

# exact score threshold for site P < 1e-4 (discretised score distribution)
scale = 100
iLO = np.round(LO * scale).astype(int)
dist = {0: 1.0}
for i in range(W):
    nd = {}
    for s, p in dist.items():
        for j in range(4):
            nd[s + iLO[i, j]] = nd.get(s + iLO[i, j], 0) + p * bg[j]
    dist = nd
sc = np.array(sorted(dist)); pr = np.array([dist[s] for s in sc])
tail = np.cumsum(pr[::-1])[::-1]
THR = sc[np.argmax(tail < P_SITE)]                          # smallest integer score with P(score>=THR) < 1e-4
idx = {b: i for i, b in enumerate(BASES)}


def encode(seqlist):
    L = max(len(q) for q in seqlist)
    X = np.full((len(seqlist), L), -1, dtype=np.int8)
    for r, q in enumerate(seqlist):
        X[r, :len(q)] = [idx.get(c, -1) for c in q]
    return X


def hit_matrix(X):
    """Boolean (n_seq x n_pos) matrix of sense-strand site starts (score >= THR)."""
    n, L = X.shape
    npos = L - W + 1
    S = np.zeros((n, npos), dtype=np.int32)
    bad = np.zeros((n, npos), dtype=bool)
    for i in range(W):
        col = X[:, i:i + npos]
        bad |= col < 0
        S += np.where(col >= 0, iLO[i][np.clip(col, 0, 3)], 0)
    return (S >= THR) & ~bad


def sites(s):
    return list(np.nonzero(hit_matrix(encode([s]))[0])[0])


def markov_shuffle_hits(seqlist, n):
    """Fraction of first-order Markov surrogates (n per sequence) containing >= 1 site."""
    hits = tot = 0
    for q in seqlist:
        q = "".join(c for c in q if c in idx)
        x = np.array([idx[c] for c in q])
        T = np.ones((4, 4)); np.add.at(T, (x[:-1], x[1:]), 1); T /= T.sum(1, keepdims=True)
        C = np.cumsum(T, 1)
        Y = np.empty((n, len(x)), dtype=np.int8); Y[:, 0] = x[0]
        U = RNG.random((n, len(x)))
        for j in range(1, len(x)):
            Y[:, j] = (U[:, j, None] > C[Y[:, j - 1]]).sum(1)
        hits += hit_matrix(Y).any(1).sum(); tot += n
    return hits / tot


# ---- set membership --------------------------------------------------------
m = master_table()
utr = pd.read_excel(TE_FILE, sheet_name="root-ctr-0h (2)", usecols=["#geneID", "5utrLen"]).rename(
    columns={"#geneID": "gene_model", "5utrLen": "utr5_len"})
m = m.merge(utr, on="gene_model", how="left").set_index("gene_model")
both = lambda t: (m.log2_te_root > t) & (m.log2_te_shoot > t)
te10 = set(m.index[both(1.0) & (m.seq_set != "")]) - set(seqs["high"])
te05 = set(m.index[both(0.5) & (m.seq_set != "")]) - set(seqs["high"])
lm = seqs["lowmid"]
sets = {
    "High-TE (>1.5; discovery set)": seqs["high"],
    "1.0 < TE ≤ 1.5 (independent)": {k: lm[k] for k in sorted(te10) if k in lm},
    "0.5 < TE ≤ 1.5 (independent)": {k: lm[k] for k in sorted(te05) if k in lm},
    "Low/mid-TE (remaining)": {k: v for k, v in lm.items() if k not in te05},
    "Low-TE (< −1.5)": seqs["low"],
}

rows, site_rows = [], []
for name, d in sets.items():
    genes = list(d); H = hit_matrix(encode([d[g] for g in genes]))
    n_hit = int(H.any(1).sum())
    n_rx = sum(bool(re.search("CTATA[AT]A[AT]A", d[g][:ATG])) for g in genes)
    for r, g in enumerate(genes):
        ulen = m.utr5_len.get(g, np.nan)
        for p in np.nonzero(H[r])[0]:
            centre = p + W // 2 - ATG
            where = "5'UTR" if (not np.isnan(ulen) and p >= ATG - ulen and p + W <= ATG) else (
                "CDS" if p >= ATG else "upstream of annotated TSS")
            site_rows.append(dict(set=name, gene_model=g, start=int(p) - ATG, centre=int(centre), utr5_len=ulen, location=where))
    nshuf = 20 if len(d) < 2000 else 3
    shuf = markov_shuffle_hits([d[g] for g in genes], nshuf)
    rows.append(dict(set=name, n=len(d), with_site=n_hit, pct=100 * n_hit / len(d), regex_exact=n_rx,
                     regex_pct=100 * n_rx / len(d), shuffled_pct=100 * shuf,
                     median_5utr_len=float(np.nanmedian([m.utr5_len.get(g, np.nan) for g in d]))))
res = pd.DataFrame(rows)
ref_low = res.iloc[4]; ref_lm = res.iloc[3]
for i, r in res.iterrows():
    for lab, ref in (("vs_low", ref_low), ("vs_lowmid", ref_lm)):
        if r.set == ref.set:
            res.loc[i, f"OR_{lab}"] = np.nan; res.loc[i, f"P_{lab}"] = np.nan; continue
        f = fisher_exact([[r.with_site, r.n - r.with_site], [ref.with_site, ref.n - ref.with_site]])
        res.loc[i, f"OR_{lab}"] = f.statistic; res.loc[i, f"P_{lab}"] = f.pvalue
    f = fisher_exact([[r.with_site, r.n - r.with_site], [round(r.shuffled_pct * r.n / 100), r.n - round(r.shuffled_pct * r.n / 100)]])
    res.loc[i, "P_vs_shuffled"] = f.pvalue
sites_df = pd.DataFrame(site_rows)
loc = sites_df.groupby("set").location.value_counts().unstack(fill_value=0)
res = res.merge(loc, left_on="set", right_index=True, how="left")
RESULTS.mkdir(exist_ok=True)
res.to_csv(RESULTS / "motif_enrichment.tsv", sep="\t", index=False)
sites_df.to_csv(RESULTS / "motif_sites.tsv", sep="\t", index=False)
with open(RESULTS / "motif_enrichment.tsv", "a") as fh:
    fh.write(f"# PWM threshold: log-odds score >= {THR/scale:.2f} (site P < {P_SITE}); background {dict(zip(BASES, bg.round(3)))}; "
             f"MEME discovery E-value of the motif = 5.5e+04\n")
pd.set_option("display.width", 250)
print(res.round(4).to_string())
print("\nsite centres (median, IQR) by set:")
print(sites_df.groupby("set").centre.describe()[["count", "25%", "50%", "75%"]])

# ---- primed genes: -150..-1 window vs the same window of the other sets ------
primed = read_fasta(MOTIF / "primed_701_upstream150.fasta")
up_sets = {"Primed (701 nuclear)": primed, **{k: {g: s[:ATG] for g, s in v.items()} for k, v in sets.items()}}
utr_g = m.utr5_len.copy(); utr_g.index = utr_g.index.str.split(".").str[0]; utr_g = utr_g[~utr_g.index.duplicated()]
prow = []
for name, d in up_sets.items():
    genes = list(d); H = hit_matrix(encode([d[g] for g in genes]))
    ul = np.array([utr_g.get(g.split(".")[0], np.nan) for g in genes])
    inutr = up = 0
    for r in range(len(genes)):
        for p in np.nonzero(H[r])[0]:
            if not np.isnan(ul[r]) and p >= ATG - ul[r]: inutr += 1
            else: up += 1
    prow.append(dict(set=name, n=len(d), with_site=int(H.any(1).sum()), pct=100 * H.any(1).mean(),
                     regex=sum(bool(re.search("CTATA[AT]A[AT]A", d[g])) for g in genes),
                     sites_in_5utr=inutr, sites_upstream_TSS=up, pct_sites_upstream_TSS=100 * up / max(inutr + up, 1),
                     median_5utr_len=float(np.nanmedian(ul)), pct_no_5utr=100 * np.mean(np.nan_to_num(ul, nan=0) == 0)))
pres = pd.DataFrame(prow)
pres["regex_pct"] = 100 * pres.regex / pres.n
for ref in ("Low/mid-TE (remaining)", "Low-TE (< −1.5)"):
    r0 = pres[pres.set == ref].iloc[0]; tag = "lowmid" if ref.startswith("Low/") else "low"
    for i, r in pres.iterrows():
        f = fisher_exact([[r.with_site, r.n - r.with_site], [r0.with_site, r0.n - r0.with_site]])
        pres.loc[i, f"OR_vs_{tag}"] = f.statistic; pres.loc[i, f"P_vs_{tag}"] = f.pvalue
pres.to_csv(RESULTS / "motif_primed_upstream150.tsv", sep="\t", index=False)
print(pres.round(4).to_string())

# ---- supplementary figure S1 (pastel) ----------------------------------------
if __name__ == "__main__":
    plt.rcParams.update({"font.family": "Liberation Sans", "font.size": 7, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.linewidth": 0.6, "pdf.fonttype": 42})
    CORAL, LAV, SAGE, PERI, BEIGE, SLATE, INK, INK2, PLUM = ("#e4876f", "#9d86cf", "#6fb39a", "#8fa8dc", "#ddd6ca",
                                                            "#b3bccd", "#0b0b0b", "#52514e", "#5b4a66")
    cols = [CORAL, LAV, "#c4b5e6", BEIGE, PERI]
    labs = ["High-TE", "1.0–1.5", "0.5–1.5", "Low/mid", "Low-TE"]
    fig = plt.figure(figsize=(7.09, 5.0))
    # panel order follows the order of citation in the text: a prevalence, b primed genes, c positions, d location
    ax_a = fig.add_axes([0.075, 0.62, 0.40, 0.31]); ax_d = fig.add_axes([0.585, 0.62, 0.40, 0.31])
    ax_b = fig.add_axes([0.09, 0.13, 0.385, 0.31]); ax_c = fig.add_axes([0.66, 0.13, 0.325, 0.31])
    a = ax_a; x = np.arange(len(res))
    a.bar(x, res.pct, color=cols, width=0.62, lw=0)
    a.scatter(x, res.shuffled_pct, marker="_", s=120, color=INK, zorder=3, label="Shuffled control")
    for i, r in res.iterrows():
        a.text(i, r.pct + 1.2, f"{r.pct:.1f}%\n({r.with_site}/{r.n})", ha="center", fontsize=5.6, color=INK)
    a.set_xticks(x); a.set_xticklabels(["High-TE\n(discovery)", "1.0–1.5\n(indep.)", "0.5–1.5\n(indep.)", "Low/mid", "Low-TE"], fontsize=6)
    a.set_ylabel("Sequences with site (%)"); a.set_ylim(0, max(res.pct) * 1.45)
    a.legend(frameon=False, loc="upper right", fontsize=6)
    a.set_title("a  Motif prevalence, −150 to +50 nt (PWM, P < 10$^{-4}$)", loc="left", fontsize=7.5)
    b = ax_b
    for i, (name, colr) in enumerate(zip(res.set, cols)):
        v = sites_df[sites_df.set == name].centre.values
        if len(v) > 300:
            v = RNG.choice(v, 300, replace=False)
        b.scatter(v, i + RNG.uniform(-0.22, 0.22, len(v)), s=5, color=colr if colr != BEIGE else "#b8ad9c", lw=0, alpha=0.9)
        med = sites_df[sites_df.set == name].centre.median()
        b.plot([med, med], [i - 0.32, i + 0.32], color=INK, lw=1.2)
    b.axvline(0, color=INK2, lw=0.6); b.text(3, -0.45, "AUG", fontsize=6, color=INK2)
    b.set_yticks(range(len(labs))); b.set_yticklabels(labs, fontsize=6); b.invert_yaxis()
    b.set_xlim(-155, 55); b.set_xlabel("Site centre relative to AUG (nt)")
    b.set_title("c  Position of sites (bar, median)", loc="left", fontsize=7.5)
    c = ax_c
    sub = res.set_index("set")[[k for k in ("5'UTR", "upstream of annotated TSS", "CDS") if k in res.columns]].fillna(0)
    frac = sub.div(sub.sum(axis=1).replace(0, np.nan), axis=0).fillna(0) * 100
    left = np.zeros(len(frac))
    for k, colr in zip(frac.columns, (SAGE, SLATE, BEIGE)):
        c.barh(np.arange(len(frac)), frac[k], left=left, color=colr, height=0.6, lw=0,
               label=k.replace("upstream of annotated TSS", "upstream of TSS"))
        left += frac[k].values
    c.set_yticks(np.arange(len(frac))); c.set_yticklabels(labs, fontsize=6)
    c.invert_yaxis(); c.set_xlabel("Sites (%)"); c.set_xlim(0, 100)
    c.legend(frameon=False, fontsize=6, loc="upper center", bbox_to_anchor=(0.45, -0.22), ncol=3)
    c.set_title("d  Site location (annotated 5′UTR)", loc="left", fontsize=7.5)
    d = ax_d
    order = ["Primed (701 nuclear)", "High-TE (>1.5; discovery set)", "0.5 < TE ≤ 1.5 (independent)",
             "Low/mid-TE (remaining)", "Low-TE (< −1.5)"]
    pr = pres.set_index("set").loc[order]
    xx = np.arange(len(pr))
    in_utr = pr.pct * pr.sites_in_5utr / (pr.sites_in_5utr + pr.sites_upstream_TSS)
    d.bar(xx, in_utr, color=SAGE, width=0.62, lw=0, label="sites in 5′UTR")
    d.bar(xx, pr.pct - in_utr, bottom=in_utr, color=SLATE, width=0.62, lw=0, label="sites upstream of TSS")
    d.bar(xx[:1], pr.pct.iloc[:1], fill=False, edgecolor=PLUM, lw=1.1, width=0.62)
    for i, (nm, r) in enumerate(pr.iterrows()):
        d.text(i, r.pct + 1.0, f"{r.pct:.1f}%\n({int(r.with_site)}/{int(r.n)})", ha="center", fontsize=5.6, color=INK)
    d.set_xticks(xx)
    d.set_xticklabels([f"{l}\n5′UTR {v:.0f} nt" for l, v in zip(["Primed (701)", "High-TE", "0.5–1.5", "Low/mid", "Low-TE"],
                                                              pr.median_5utr_len)], fontsize=6)
    d.set_ylabel("Sequences with site (%)"); d.set_ylim(0, max(pr.pct) * 1.35); d.set_xlim(-0.6, len(pr) - 0.4)
    d.legend(frameon=False, fontsize=6, loc="upper right", bbox_to_anchor=(1.0, 1.02))
    d.set_title("b  Translationally primed genes, −150 to −1 nt", loc="left", fontsize=7.5)
    for ext in ("png", "pdf"):
        fig.savefig(FIGURES / f"FigureS1_motif_reanalysis.{ext}", dpi=600)
    print("figure written")
