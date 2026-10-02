# Translational priming in *Arabidopsis thaliana*

Code and derived data to reproduce the analyses, figures and Supporting Information tables of
*"Translational priming …"* (Vélez-Bermúdez *et al.*, in preparation, *New Phytologist*).

The analysis integrates matched RNA-seq and Ribo-seq data from roots and shoots
(Vélez-Bermúdez *et al.*, 2025, *Scientific Data* 12: 504) with the root and shoot
proteomes of Jain & Schmidt (2023, *Molecular & Cellular Proteomics* 23: 100685).

## Definitions (all in `scripts/te_common.py`)

| Term | Definition |
|---|---|
| TE | Ribo-seq RPKM / RNA-seq RPKM, per tissue. Needs ≥ 1 RNA-seq read (undefined if RNA-seq RPKM = 0); TE = 0 if there are no Ribo-seq reads. No higher read-count threshold is applied |
| high-TE / mid-TE / low-TE | log2TE > 1.5 / −1.5 ≤ log2TE ≤ 1.5 / log2TE < −1.5 (TE = 0 counted as low-TE and flagged) |
| Translationally primed | log2 RNA-seq RPKM < 0 **and** log2TE > 1.5 in the same tissue (hence ≥ 1 RNA-seq and ≥ 1 Ribo-seq read) |
| Canonical high-TE | log2TE > 1.5 with log2 RNA-seq RPKM ≥ 0 in every high-TE tissue |
| Sequence-analysis sets | log2 RNA-seq RPKM ≥ 2 in both tissues (13,133 genes); high-TE (27) / low-TE (1,335) in both tissues; remainder low/mid-TE (11,720) |
| Motif validation sets | 1.0 < log2TE ≤ 1.5 (60) and 0.5 < log2TE ≤ 1.5 (625) in both tissues, taken from the low/mid-TE set (not used for motif discovery) |
| Protein evidence | protein detected with High or Medium FDR confidence in the root or shoot proteome (any pH), matched by AGI locus |

## Input files

| File | Source |
|---|---|
| `data/raw/TranslationEfficiency2.xlsx` | Processed per-tissue RPKM values (W.-D. Lin), derived from Vélez-Bermúdez *et al.* (2025) — *not distributed; place here* |
| `data/raw/Jain_Schmidt_2023_mmc3.xlsx` | Supplementary dataset of Jain & Schmidt (2023), doi:10.1016/j.mcpro.2023.100685 — *not distributed; place here* |
| `data/raw/TE_full_annotation.xlsx` | Araport11 descriptions and functional categories — *not distributed; place here* |
| `data/raw/motif_analysis_Supplementary_Data_5.xlsx` | Original MEME analysis (motifs discovered in low-TE sequences) — *not distributed; place here* |
| `data/raw/KEGG_ath00999_genes.txt` | Genes of KEGG pathway ath00999 (from rest.kegg.jp/link/ath/path:ath00999) |
| `data/motif/*_200bp.fasta` | −150/+50 nt around the annotated AUG (A of AUG at 0-based index 150) for the high-TE, low-TE and low/mid-TE sets (TAIR10 genome, Araport11 models) |
| `data/motif/primed_701_upstream150.fasta` | −150/−1 nt upstream of the AUG of the 701 nuclear translationally primed genes |
| `data/motif/CTATAWAWA_motif.meme` | MEME motif 6 (CTATA[A/T]A[A/T]A-like), discovered in the 27 high-TE sequences |
| `data/gprofiler/gene_lists/*.txt` | Gene lists submitted to g:Profiler (regenerated and checked by script 05) |
| `data/gprofiler/results/` | g:Profiler output (version e114_eg62_p19_27110d83, run 2026-09-28) |

## Run

```bash
pip install -r requirements.txt
cd scripts
python 01_verify_manuscript_numbers.py   # master table + every number quoted in the text -> results/
python 02_make_figures.py                # Figures 1 and 2 (PNG 600 dpi + PDF) and source data -> figures/
python 04_motif_analysis.py              # motif re-analysis, incl. primed genes -> results/motif_*.tsv, Fig. S1
python 05_gprofiler_enrichment.py        # writes/checks g:Profiler gene lists (add --run to re-query g:Profiler; internet needed)
python 03_supporting_tables.py           # Tables S1–S9 and Table 1 -> supporting_information/, results/
```

Order matters: `03` uses the outputs of `01`, `02` (Fig. 2b source data), `04` (motif tables) and the
g:Profiler results in `data/gprofiler/results/`.

## Methods summary

* **Protein evidence vs read support** – primed genes vs RNA-low genes that are not primed, at matched minimum
  Ribo-seq read support (≥ 1, 5, 10, 20 reads) per tissue; two-sided Fisher's exact tests (Table 1, Fig. 1e,f).
* **Enrichment** – g:Profiler g:GOSt; GO:BP, GO:MF, GO:CC (electronic annotations excluded) and KEGG;
  Benjamini–Hochberg FDR < 0.05; annotated genome as background; the primed set is also run against the 24,755
  nuclear genes detected by RNA-seq. ath00999 is additionally tested with a one-sided hypergeometric test.
* **Motif** – sense-strand PWM scan (pseudocount 0.1; background = composition of all windows); sites with
  P < 1e-4 from the exact score distribution (as in FIMO); Fisher's exact tests vs remaining low/mid-TE and
  low-TE genes; first-order Markov (dinucleotide) shuffles as control (20 per sequence; 3 for sets > 2,000);
  sites located relative to the annotated 5′UTR (Araport11 5′UTR length). Statistics for the discovery set are
  circular and are reported for reference only.

## Outputs

* `results/master_gene_table.tsv` – one row per gene model with all derived values
* `results/manuscript_numbers.tsv` – every quantity cited in the manuscript
* `results/Table1_protein_evidence.tsv`
* `results/motif_enrichment.tsv`, `results/motif_primed_upstream150.tsv`, `results/motif_sites.tsv`
* `figures/Figure1_*`, `figures/Figure2_*`, `figures/FigureS1_*` (PNG + PDF) and `Fig*.tsv` source data
* `supporting_information/Supporting_Information_Tables_S1-S9.xlsx`

## Software

Python 3.11; pandas 3.0.2; NumPy 2.4.4; SciPy 1.17.1; Matplotlib 3.10.9; openpyxl 3.1.5;
gprofiler-official (g:Profiler version e114_eg62_p19_27110d83).
De novo motif discovery (MEME, MEME Suite) was run separately; its output is not regenerated here.

## License

Code: MIT License. Derived data tables: CC BY 4.0.
