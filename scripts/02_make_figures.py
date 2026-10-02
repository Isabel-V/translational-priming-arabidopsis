"""Composite Figures 1 and 2 (proposed) for the translational-priming Letter.

All values are computed from the raw tables (see te_common.py), except the
functional-category summary (Fig. 1e), which is read from the functional
annotation workbook used for the manuscript text (TE_full_annotation.xlsx,
sheet Summary_ByCategory), and the start-codon nucleotide matrices (Fig. 2c),
which are read from Wen-Dar Lin's ATG_region_matrixes sheet.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy.stats import fisher_exact
from te_common import master_table, FIGURES, ANNOT_FILE, TE_FILE, TE_HIGH, TE_LOW

FIGURES.mkdir(exist_ok=True)
plt.rcParams.update({
    "font.family": "Liberation Sans", "font.size": 7, "axes.titlesize": 7.5, "axes.labelsize": 7,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.2, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42, "svg.fonttype": "none",
    # math text in the same typeface, so that $n$ and gene identifiers are set in Liberation Sans italic
    "mathtext.fontset": "custom", "mathtext.rm": "Liberation Sans", "mathtext.it": "Liberation Sans:italic",
    "mathtext.bf": "Liberation Sans:bold", "mathtext.default": "regular",
})


def fmt_agi(a):
    """Arabidopsis locus identifier in journal style: At2g31141 (lower case except A)."""
    a = a.split(".")[0]
    return "At" + a[2:].lower() if a.upper().startswith("AT") else a


def it(s):
    """Italic text via mathtext (keeps Liberation Sans)."""
    return r"$\it{" + s.replace("'", r"\prime") + "}$"
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#898781", "#e1e0d9"
HIGH, LOW, MID = "#e4876f", "#8fa8dc", "#ddd6ca"      # pastel diverging: coral (high/primed), periwinkle (low), warm beige (mid)
CTRL = "#b3bccd"                                      # reference / control group
S1, S2 = "#9d86cf", "#6fb39a"                         # categorical slots 1-2 (Kozak A-3 / G+4)
PRIMED_FILL = "#fdf0eb"
PEACH, LAV, PLUM, CORAL_TXT, BEIGE_BAR, BAND = "#f5c6b8", "#9d86cf", "#5b4a66", "#b0553f", "#d9d2c7", "#f8f5f0"

m = master_table()


def panel_label(ax, s, x=None, y=None):
    b = ax.get_position()
    fx = b.x0 - 0.065 if x is None else x
    ax.figure.text(fx, b.y1 + 0.012, s, fontsize=9, fontweight="bold", va="bottom", ha="left", color=INK)


def scatter(ax, t, labels):
    d = m[m[f"rna_detected_{t}"] & m[f"ribo_detected_{t}"]].copy()
    x, y = d[f"log2_rna_{t}"], d[f"log2_te_{t}"]
    ax.add_patch(plt.Rectangle((-11, TE_HIGH), 11, 17, color=PRIMED_FILL, zorder=0, lw=0))
    cls = d[f"te_class_{t}"]
    for c, col, s, a in (("mid-TE", MID, 1.2, 0.55), ("low-TE", LOW, 1.2, 0.5), ("high-TE", PEACH, 2.4, 0.95)):
        k = (cls == c) & ~d[f"primed_{t}"]
        ax.scatter(x[k], y[k], s=s, c=col, lw=0, alpha=a, rasterized=True, zorder=1)
    k = d[f"primed_{t}"]
    ax.scatter(x[k], y[k], s=5, c=HIGH, edgecolors=PLUM, linewidths=0.3, zorder=3, rasterized=True)
    ax.axhline(TE_HIGH, color=INK2, lw=0.5, ls=(0, (3, 2)), zorder=2)
    ax.axhline(TE_LOW, color=INK2, lw=0.5, ls=(0, (3, 2)), zorder=2)
    ax.axvline(0, color=INK2, lw=0.5, ls=(0, (1, 1.5)), zorder=2)
    ax.set_xlim(-10.5, 12.5); ax.set_ylim(-11, 18)
    ax.set_xlabel(r"RNA-seq, log$_2$ RPKM"); ax.set_ylabel(r"Translational efficiency, log$_2$ TE")
    n_primed = int(k.sum())
    ax.text(-10.2, 17.7, f"Translationally primed ($\\it{{n}}$ = {n_primed})", color=CORAL_TXT, fontsize=6.3, va="top", fontweight="bold")
    ax.set_title(f"{t.capitalize()}s ($\\it{{n}}$ = {len(d):,} genes)", loc="left", color=INK)
    for name, agi, dx, dy in labels:
        r = d[d.agi == agi]
        if r.empty:
            continue
        r = r.iloc[0]
        ax.annotate(fmt_agi(name), (r[f"log2_rna_{t}"], r[f"log2_te_{t}"]), xytext=(r[f"log2_rna_{t}"] + dx, r[f"log2_te_{t}"] + dy),
                    fontsize=5.8, color=INK, fontstyle="italic",
                    arrowprops=dict(arrowstyle="-", lw=0.4, color=INK2, shrinkA=0, shrinkB=1.5), zorder=4)


# =============================== FIGURE 1 ===============================
fig = plt.figure(figsize=(7.09, 6.3))
# panel order follows the order of citation in the text: a gene sets, b/c scatter, d categories, e/f protein evidence
axC = fig.add_axes([0.165, 0.575, 0.155, 0.37]); axA = fig.add_axes([0.425, 0.575, 0.235, 0.37]); axB = fig.add_axes([0.745, 0.575, 0.235, 0.37])
axE = fig.add_axes([0.205, 0.085, 0.115, 0.37]); axD1 = fig.add_axes([0.425, 0.085, 0.235, 0.37]); axD2 = fig.add_axes([0.745, 0.085, 0.235, 0.37])

scatter(axA, "root", [("S8H", "AT3G12900", 2.2, 2.2), ("CYP82C4", "AT4G31940", 2.4, 0.2), ("CRK28", "AT4G21400", 1.8, -1.6),
                      ("AT3G55860", "AT3G55860", 2.0, -2.2)])
scatter(axB, "shoot", [("IRT1", "AT4G19690", -2.5, 3.0), ("AT2G31141", "AT2G31141", 1.6, 1.4), ("S8H", "AT3G12900", 3.5, -4.2)])
panel_label(axA, "b"); panel_label(axB, "c")
leg = [Line2D([], [], marker="o", ls="", ms=3, mfc=HIGH, mec=PLUM, mew=0.3, label="Primed (log$_2$ mRNA < 0, log$_2$TE > 1.5)"),
       Line2D([], [], marker="o", ls="", ms=3, mfc=PEACH, mec="none", label="Canonical high-TE"),
       Line2D([], [], marker="o", ls="", ms=3, mfc=MID, mec="none", label="Mid-TE"),
       Line2D([], [], marker="o", ls="", ms=3, mfc=LOW, mec="none", label="Low-TE")]
axB.legend(handles=leg, loc="lower right", frameon=False, handletextpad=0.2, borderaxespad=0.1, fontsize=5.8)

# a: gene-set cascade
det = (m.rna_detected_root | m.rna_detected_shoot) & (m.ribo_detected_root | m.ribo_detected_shoot)
split = m.primed_any & (((m.te_class_root == "high-TE") & ~m.primed_root) | ((m.te_class_shoot == "high-TE") & ~m.primed_shoot))
rows = [("Gene models", len(m), BEIGE_BAR), ("Detected (RNA + Ribo)", int(det.sum()), BEIGE_BAR),
        ("High-TE (≥ 1 tissue)", int(m.high_te_any.sum()), PEACH),
        ("Primed (≥ 1 tissue)", int(m.primed_any.sum()), HIGH),
        ("  roots", int(m.primed_root.sum()), HIGH), ("  shoots", int(m.primed_shoot.sum()), HIGH),
        ("  both tissues", int(m.primed_both.sum()), HIGH),
        ("  both tissues, nuclear", int((m.primed_both & (m.organelle == "nuclear")).sum()), HIGH)]
ypos = np.arange(len(rows))[::-1]
for yy, (lab, v, col) in zip(ypos, rows):
    axC.barh(yy, v, color=col, height=0.62, lw=0)
    axC.text(v * 1.12, yy, f"{v:,}", va="center", fontsize=6.2, color=INK)
axC.set_yticks(ypos); axC.set_yticklabels([r[0] for r in rows], fontsize=6.2)
axC.set_xscale("log"); axC.set_xlim(20, 2.2e5); axC.set_xlabel("Number of genes (log scale)")
axC.tick_params(axis="y", length=0)
pct = 100 * m.primed_any.sum() / m.high_te_any.sum()

axC.set_title("Gene sets", loc="left", color=INK)
panel_label(axC, "a", x=0.008)

# e, f: protein-level evidence vs Ribo-seq read support, RNA-low genes only
thr = [1, 5, 10, 20]
stats_rows = []
for ax, t in ((axD1, "root"), (axD2, "shoot")):
    base = m[m[f"ribo_detected_{t}"] & (m[f"log2_rna_{t}"] < 0)]
    yp, yc, ns = [], [], []
    for k in thr:
        a = base[base[f"primed_{t}"] & (base[f"ribo_reads_{t}"] >= k)]
        b = base[~base[f"primed_{t}"] & (base[f"ribo_reads_{t}"] >= k)]
        f = fisher_exact([[a.protein_any.sum(), len(a) - a.protein_any.sum()], [b.protein_any.sum(), len(b) - b.protein_any.sum()]])
        yp.append(100 * a.protein_any.mean()); yc.append(100 * b.protein_any.mean()); ns.append((len(a), len(b), f.statistic, f.pvalue))
        stats_rows.append(dict(tissue=t, min_ribo_reads=k, primed_n=len(a), primed_protein=int(a.protein_any.sum()),
                               primed_pct=round(100 * a.protein_any.mean(), 1), control_n=len(b),
                               control_protein=int(b.protein_any.sum()), control_pct=round(100 * b.protein_any.mean(), 1),
                               odds_ratio=round(f.statistic, 2), fisher_p=f.pvalue))
    xs = np.arange(len(thr)); w = 0.36
    ax.bar(xs - w / 2 - 0.02, yp, w, color=HIGH, lw=0, label="Primed")
    ax.bar(xs + w / 2 + 0.02, yc, w, color=CTRL, lw=0, label="RNA-low, not primed")
    for i, (npri, nctl, orr, p) in enumerate(ns):
        star = "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "n.s."
        ax.text(xs[i], max(yp[i], yc[i]) + 2.5, star, ha="center", fontsize=6.2, color=INK)
        ax.text(xs[i] - w / 2 - 0.02, 1.5, f"{npri}", ha="center", va="bottom", fontsize=5, color=INK2, rotation=90)
        ax.text(xs[i] + w / 2 + 0.02, 1.5, f"{nctl}", ha="center", va="bottom", fontsize=5, color=INK2, rotation=90)
    ax.set_xticks(xs); ax.set_xticklabels([f"≥ {k}" for k in thr]); ax.set_ylim(0, 68)
    ax.set_xlabel("Minimum Ribo-seq reads"); ax.set_ylabel("Genes with protein evidence (%)")
    ax.yaxis.grid(True, color=GRID, lw=0.4); ax.set_axisbelow(True)
    ax.set_title(f"{t.capitalize()}s, genes with log$_2$ mRNA < 0", loc="left", color=INK)
axD1.legend(frameon=False, loc="upper left", handlelength=1, borderaxespad=0.1)
axD2.legend(frameon=False, loc="upper left", handlelength=1, borderaxespad=0.1)
panel_label(axD1, "e"); panel_label(axD2, "f")
pd.DataFrame(stats_rows).to_csv(FIGURES / "Fig1ef_protein_vs_read_support.tsv", sep="\t", index=False)

# d: fraction high-TE per functional category (values as cited in the text)
cat = pd.read_excel(ANNOT_FILE, sheet_name="Summary_ByCategory", header=2).dropna(subset=["Category"])
cat["high"] = cat["High-TE %"].astype(str).str.rstrip("%").astype(float)
ren = {"Energy metabolism (Plastid-encoded)": "Energy, plastid-enc.",
       "Energy metabolism (Mitochondrial-encoded)": "Energy, mito-enc.",
       "Energy metabolism (Nuclear-encoded)": "Energy, nuclear-enc.",
       "Ribosomal/translation machinery": "Ribosome/translation", "Uncharacterized/unknown function": "Uncharacterized",
       "Other characterized": "Other characterized", "Transcription factor": "Transcription factor"}
cat["lab"] = cat.Category.map(lambda s: ren.get(s, s))
cat = cat.sort_values("high")
axE.barh(np.arange(len(cat)), cat.high, color=[HIGH if v > 10 else PEACH for v in cat.high], height=0.62, lw=0)
for i, (v, n) in enumerate(zip(cat.high, cat["Total genes"])):
    axE.text(v + 2, i, f"{v:.1f}%", va="center", fontsize=5.6, color=INK)
axE.axvline(3.4, color=INK2, lw=0.6, ls=(0, (3, 2)))
axE.text(4.2, len(cat) - 0.35, "genome-wide 3.4%", fontsize=5.6, color=INK2, va="bottom")
axE.set_yticks(np.arange(len(cat))); axE.set_yticklabels([f"{l} ({int(n):,})" for l, n in zip(cat.lab, cat["Total genes"])], fontsize=5.8); axE.tick_params(axis="y", length=0)
axE.set_xlim(0, 100); axE.set_xlabel("High-TE genes (%)")
axE.set_title("Functional categories", loc="left", color=INK)
panel_label(axE, "d", x=0.008)

for ext in ("png", "pdf"):
    fig.savefig(FIGURES / f"Figure1_translational_priming_landscape.{ext}", dpi=600)
plt.close(fig)

# =============================== FIGURE 2 ===============================
iron = [  # (label, AGI, group)
    ("IRT1", "AT4G19690", "Transport"), ("IRT2", "AT4G19680", "Transport"), ("FRO2", "AT1G01580", "Transport"),
    ("MTP8", "AT3G58060", "Transport"), ("VTL3", "AT3G43630", "Transport"), ("VTL2", "AT1G76800", "Transport"),
    ("S8H", "AT3G12900", "Coumarin & nicotianamine synthesis"), ("CYP82C4", "AT4G31940", "Coumarin & nicotianamine synthesis"),
    ("CYP82C3", "AT4G31950", "Coumarin & nicotianamine synthesis"), ("F6'H1", "AT3G13610", "Coumarin & nicotianamine synthesis"),
    ("NAS2", "AT5G56080", "Coumarin & nicotianamine synthesis"), ("NAS3", "AT1G09240", "Coumarin & nicotianamine synthesis"),
    ("FIT", "AT2G28160", "Transcriptional regulators"), ("bHLH38", "AT3G56970", "Transcriptional regulators"),
    ("bHLH39", "AT3G56980", "Transcriptional regulators"), ("bHLH100", "AT2G41240", "Transcriptional regulators"),
    ("bHLH101", "AT5G04150", "Transcriptional regulators"), ("PYE", "AT3G47640", "Transcriptional regulators"),
    ("MYB72", "AT1G56160", "Transcriptional regulators"), ("BTS", "AT3G18290", "Transcriptional regulators"),
]
fig = plt.figure(figsize=(7.09, 5.0))
# panel order follows the order of citation in the text: a small ORFs, b iron network, c start-codon context
gs = fig.add_gridspec(2, 2, width_ratios=[1, 1.25], height_ratios=[1, 1], hspace=0.62, wspace=0.34,
                      left=0.075, right=0.975, top=0.93, bottom=0.13)
axA = fig.add_subplot(gs[:, 1]); axB = fig.add_subplot(gs[0, 0]); axC = fig.add_subplot(gs[1, 0])

rows = []
yy, ylab, yv, groups = 0, [], [], []
for g in ("Transport", "Coumarin & nicotianamine synthesis", "Transcriptional regulators"):
    items = [i for i in iron if i[2] == g]
    groups.append((g, yy, yy + len(items) - 1))
    for lab, agi, _ in items:
        r = m[m.agi == agi].iloc[0]
        for t, mk, off in (("root", "o", -0.17), ("shoot", "s", 0.17)):
            v = r[f"log2_te_{t}"]
            if pd.isna(v):
                continue
            primed = bool(r[f"primed_{t}"]); reads = r[f"ribo_reads_{t}"]
            col = HIGH if primed else (MID if v <= TE_HIGH else PEACH)
            fc = col if (reads >= 10 or not primed) else "white"
            axA.scatter(v, yy + off, marker=mk, s=16, facecolors=fc, edgecolors=col if not primed else PLUM, linewidths=0.6, zorder=3)
            rows.append(dict(gene=lab, agi=agi, group=g, tissue=t, log2_rna=round(r[f"log2_rna_{t}"], 2), log2_te=round(v, 2),
                             primed=primed, ribo_reads=int(reads), rna_reads=int(r[f"rna_reads_{t}"]), protein=bool(r.protein_any)))
        ylab.append(lab + (" †" if r.protein_any else "")); yv.append(yy); yy += 1
    yy += 0.8
axA.axvline(TE_HIGH, color=INK2, lw=0.5, ls=(0, (3, 2))); axA.axvline(0, color=GRID, lw=0.6)
axA.set_yticks(yv); axA.set_yticklabels(ylab, fontsize=6.2); axA.invert_yaxis(); axA.tick_params(axis="y", length=0)
for lab in axA.get_yticklabels():
    lab.set_fontstyle("italic")
axA.set_xlabel(r"Translational efficiency, log$_2$ TE"); axA.set_xlim(-1.8, 7.2)
for g, a, b in groups:
    axA.text(7.1, a - 0.75, g, fontsize=6.4, fontweight="bold", color=INK, ha="right")
    axA.axhspan(a - 0.5, b + 0.5, color=BAND if g != "Transcriptional regulators" else "white", zorder=0, lw=0)
axA.set_title("Iron-acquisition network", loc="left", color=INK)
leg = [Line2D([], [], marker="o", ls="", ms=4, mfc=HIGH, mec=PLUM, mew=0.6, label="Roots, primed"),
       Line2D([], [], marker="s", ls="", ms=4, mfc=HIGH, mec=PLUM, mew=0.6, label="Shoots, primed"),
       Line2D([], [], marker="o", ls="", ms=4, mfc="white", mec=PLUM, mew=0.6, label="Primed, <10 Ribo-seq reads"),
       Line2D([], [], marker="o", ls="", ms=4, mfc=PEACH, mec=PEACH, label="Canonical high-TE"),
       Line2D([], [], marker="o", ls="", ms=4, mfc=MID, mec=MID, label="Mid-TE")]
axA.legend(handles=leg, loc="lower right", frameon=True, facecolor="white", edgecolor="none", fontsize=5.6, handletextpad=0.2,
           borderaxespad=0.2)
panel_label(axA, "b", x=axA.get_position().x0 - 0.11)
pd.DataFrame(rows).to_csv(FIGURES / "Fig2b_iron_genes.tsv", sep="\t", index=False)

# b: four shared nuclear primed genes with the highest TE (small ORFs)
s41 = m[m.primed_both & (m.organelle == "nuclear")].copy()
s41["mean_te"] = s41[["log2_te_root", "log2_te_shoot"]].mean(axis=1)
top = s41.sort_values("mean_te", ascending=False).head(4)
xs = np.arange(4); w = 0.36
axB.bar(xs - w / 2 - 0.02, top.log2_te_root, w, color=HIGH, lw=0, label="Roots")
axB.bar(xs + w / 2 + 0.02, top.log2_te_shoot, w, color=LAV, lw=0, label="Shoots")
for i, (_, r) in enumerate(top.iterrows()):
    axB.text(i - w / 2 - 0.02, r.log2_te_root + 0.3, f"{int(r.ribo_reads_root):,}", ha="center", fontsize=5, color=INK2, rotation=90, va="bottom")
    axB.text(i + w / 2 + 0.02, r.log2_te_shoot + 0.3, f"{int(r.ribo_reads_shoot):,}", ha="center", fontsize=5, color=INK2, rotation=90, va="bottom")
axB.axhline(TE_HIGH, color=INK2, lw=0.5, ls=(0, (3, 2)))
axB.set_xticks(xs); axB.set_xticklabels([fmt_agi(a) for a in top.gene_model], fontsize=5.8, fontstyle="italic")
for i, l in enumerate(top.predicted_protein_length_aa):
    axB.text(i, -0.115, f"{int(l)} aa", transform=axB.get_xaxis_transform(), ha="center", va="top", fontsize=5.8, color=INK)
axB.set_ylabel(r"log$_2$ TE"); axB.set_ylim(0, 18.5)
axB.legend(frameon=False, loc="upper right", ncol=2, handlelength=1, borderaxespad=0.1)
axB.set_title("Top-TE shared primed genes (small ORFs)", loc="left", color=INK)
panel_label(axB, "a", x=0.008)

# c: Kozak-like context across TE classes (Wen-Dar Lin, ATG_region_matrixes; columns = -3,-2,-1,A,T,G,+4,+5,+6)
mx = pd.read_excel(TE_FILE, sheet_name="ATG_region_matrixes", header=None)
blocks = []
for i, v in mx[0].items():
    if isinstance(v, str) and ("TE" in v):
        n = int(v.split("n=")[1].rstrip(")"))
        sub = mx.iloc[i + 1:i + 5].set_index(0).astype(float)
        blocks.append((v, n, sub))
order = [("TE lower than -1.5", "Low-TE\n(< −1.5)"), ("TE midLow", "Low/mid-TE"), ("TE higher than 0.5", "> 0.5"),
         ("TE higher than 1.0", "> 1.0"), ("TE higher than 1.5", "High-TE\n(> 1.5)")]
xa, ya3, gp4, ns = [], [], [], []
for key, lab in order:
    v, n, sub = [b for b in blocks if b[0].startswith(key)][0]
    xa.append(lab); ns.append(n); ya3.append(100 * sub.loc["A", 1]); gp4.append(100 * sub.loc["G", 7])
xs = np.arange(len(order))
axC.plot(xs, ya3, color=S1, lw=2, marker="o", ms=4, label="A at −3")
axC.plot(xs, gp4, color=S2, lw=2, marker="s", ms=4, label="G at +4")
axC.text(xs[-1] + 0.12, ya3[-1], f"{ya3[-1]:.0f}%", color=INK, fontsize=6, va="center")
axC.text(xs[0] - 0.12, ya3[0], f"{ya3[0]:.0f}%", color=INK, fontsize=6, va="center", ha="right")
axC.text(xs[-1] + 0.12, gp4[-1] - 3, f"{gp4[-1]:.0f}%", color=INK, fontsize=6, va="center")
axC.text(xs[0] - 0.12, gp4[0], f"{gp4[0]:.0f}%", color=INK, fontsize=6, va="center", ha="right")
axC.axhline(25, color=MUTED, lw=0.5, ls=(0, (1, 1.5))); axC.text(4.35, 26, "random", fontsize=5.4, color=MUTED, ha="right")
axC.set_xticks(xs); axC.set_xticklabels([f"{l}\n$\\it{{n}}$ = {n:,}" for l, n in zip(xa, ns)], fontsize=5.6)
axC.set_xlim(-0.6, 4.6); axC.set_ylim(15, 85); axC.set_ylabel("Nucleotide frequency (%)")
axC.set_xlabel(r"TE class (log$_2$ mRNA ≥ 2 in both tissues)")
axC.yaxis.grid(True, color=GRID, lw=0.4); axC.set_axisbelow(True)
axC.legend(frameon=False, loc="upper left", handlelength=1.4, borderaxespad=0.1)
axC.set_title("Start-codon context (Kozak-like)", loc="left", color=INK)
panel_label(axC, "c", x=0.008)
n_hi, n_lo = ns[-1], ns[0]
a_hi, a_lo = round(ya3[-1] / 100 * n_hi), round(ya3[0] / 100 * n_lo)
f = fisher_exact([[a_hi, n_hi - a_hi], [a_lo, n_lo - a_lo]])
pd.DataFrame(dict(te_class=[o[1].replace("\n", " ") for o in order], n=ns, A_minus3_pct=np.round(ya3, 1), G_plus4_pct=np.round(gp4, 1))).to_csv(
    FIGURES / "Fig2c_kozak_context.tsv", sep="\t", index=False)
with open(FIGURES / "Fig2c_kozak_context.tsv", "a") as fh:
    fh.write(f"# Fisher exact, A at -3, high-TE vs low-TE: OR={f.statistic:.2f}, P={f.pvalue:.2e}\n")

for ext in ("png", "pdf"):
    fig.savefig(FIGURES / f"Figure2_iron_network_small_ORFs_kozak.{ext}", dpi=600)
plt.close(fig)
print("figures written")
