# Metagenomic_GRE

**Cellular and vesicular metagenomes of the Guadalquivir River estuary**

Code for the Master's thesis **"Exploring the Bacteria Extracellular Vesicles fraction of the Guadalquivir River estuary"** (Alba Mata González, Máster Universitario en Análisis de Datos Ómicos y Biología de Sistemas, Universidad de Sevilla – Universidad Internacional de Andalucía).

The study compares the cellular fraction (CF) and a vesicle-enriched fraction (VF) of the microbiome along the salinity gradient of the Guadalquivir River estuary (SW Spain), using shotgun metagenomics. This repository holds the two scripts that turn the annotation tables into every figure, table and statistic of the manuscript.

## Contents

| File | Language | What it does |
| --- | --- | --- |
| `01_aggregate_orfs_by_ko.py` | Python 3 | Filters the SqueezeMeta ORF tables to prokaryotic coding sequences, sums read counts per KEGG orthologue (KO) and writes the tables read by the R script |
| `02_GRE_CF_VF_analysis.R` | R | Taxonomic, functional and nitrogen-cycle analyses; produces Figures 2–5, Table 1, Figure S1 and the supplementary tables |
| `README.md` | – | This file |

Run the scripts in the order of their numbers.

## Study design

Water was sampled on 11 November 2024 at ten stations. The CF was sequenced at all ten; the VF at four of them.

| Zone | Stations (CF) | Station with VF | Group used in the tests |
| --- | --- | --- | --- |
| Coastal Waters | 2 | 2 | Marine-influenced |
| Inlet | 5, 6, 7 | 6 | Marine-influenced |
| Transition | 9, 10, 11 | 10 | Low-salinity |
| Inner | 13, 14, 15 | 14 | Low-salinity |

Two sets of sample identifiers appear in the raw files. Both scripts keep them as they are; the R script maps them to station numbers in a single table (`DESIGN`).

| Station | mOTUs tables | SqueezeMeta tables |
| --- | --- | --- |
| 2 | `A1` / `V_A1` | `Z3_2` (relabelled `Z4_2` in R) / `Z4_V2` |
| 5, 6, 7 | `B1`, `B2`, `B3` / `V_B2` | `Z3_5`, `Z3_6`, `Z3_7` / `Z3_V6` |
| 9, 10, 11 | `C1`, `C2`, `C3` / `V_C2` | `Z2_9`, `Z2_10`, `Z2_11` / `Z2_V10` |
| 13, 14, 15 | `D1`, `D2`, `D3` / `V_D2` | `Z1_13`, `Z1_14`, `Z1_15` / `Z1_V14` |

## Requirements

**Python** 3.8 or later. Only the standard library is used.

**R** 4.5.0 (the version used for the thesis) with these packages:

| Needed | Optional (their blocks are skipped if missing) |
| --- | --- |
| tidyverse, readxl, phyloseq, vegan, patchwork, iNEXT | clusterProfiler (over-representation analysis), eulerr (Venn diagrams of Figure 4A), KEGGREST (only to rebuild the KEGG pathway tables), microbiome (Aitchison robustness check) |

The script loads the needed packages with `pacman::p_load()`, which installs missing CRAN packages. `phyloseq`, `clusterProfiler`, `KEGGREST` and `microbiome` come from Bioconductor and must be installed beforehand:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("phyloseq", "clusterProfiler", "KEGGREST", "microbiome"))
```

The exact package versions of a run are written to `outputs/sessionInfo.txt`.

## Input data

The input files are not in this repository. [State here where they are: accession number of the raw reads, or "available from the authors on request".]

Upstream processing, done before these scripts:

- **mOTUs v3** on the unassembled reads, with GTDB taxonomy.
- **SqueezeMeta** for assembly, ORF prediction, taxonomic assignment (DIAMOND against NCBI nr, last common ancestor) and KEGG annotation. The CF is one project with ten samples; the VF is four projects, one per station.
- **NCycDB** for nitrogen-cycling gene families.

## Step 1 – `01_aggregate_orfs_by_ko.py`

```bash
python3 01_aggregate_orfs_by_ko.py \
    --fc  annotation/13.FC_Plan_nacional.orftable.prok.tsv \
    --fv  VESICULAS/13.Z1_V14.orftable VESICULAS/13.Z2_V10.orftable \
          VESICULAS/13.Z3_V6.orftable  VESICULAS/13.Z4_V2.orftable \
    --kegg_txt     docs/keggfun2.txt \
    --out_tpm      resultado/kegg_tpm.tsv \
    --out_tax      resultado/kegg_tax.tsv \
    --out_kegg_ref resultado/kegg_ref.tsv
```

| Argument | Input or output | Content |
| --- | --- | --- |
| `--fc` | Input | SqueezeMeta ORF table of the CF project, already restricted to prokaryotic ORFs |
| `--fv` | Input | SqueezeMeta ORF tables of the VF, one file per station; the station name is read from the file name |
| `--kegg_txt` | Input | KEGG reference file of SqueezeMeta (`keggfun2.txt`): KO, gene name, function and pathways |
| `--out_tpm` | Output | One row per KO and fraction, with raw read counts (`Raw_<station>`) and scaled abundances (`TPM_<station>`) |
| `--out_tax` | Output | One row per unique combination of KO, fraction, station and lineage |
| `--out_kegg_ref` | Output | One row per KO and pathway, with the hierarchy split into levels |

What the script does:

1. Keeps coding sequences only (`Molecule = CDS`).
2. In the VF tables, keeps ORFs whose taxonomy starts with Bacteria or Archaea. The CF table is expected to be filtered in the same way beforehand.
3. Sums raw read counts of all ORFs that share a KO, per station.
4. Divides each sum by the mean ORF length of that KO and scales it by a denominator **shared by all samples of both fractions**.

Because of step 4, the `TPM_` columns are not normalised within each sample. The R script converts them to within-sample proportions where a per-sample composition is needed.

## Step 2 – `02_GRE_CF_VF_analysis.R`

1. Open the script and edit the `CFG` block at the top: output folder and paths of the input files.
2. Run it from top to bottom in a clean R session, or with `Rscript 02_GRE_CF_VF_analysis.R`.

### Input files

| `CFG` entry | File | Produced by |
| --- | --- | --- |
| `in_motus_counts` | `MetaQVIR_motus_filtered.txt` | mOTUs: counts per mOTU and sample |
| `in_motus_tax` | `MetaQVIR_motus_GTDB_taxonomy.tsv` | mOTUs: GTDB taxonomy |
| `in_motus_meta` | `analisis_metada.txt` | Sample metadata for the mOTUs tables |
| `in_kegg_tpm` | `kegg_tpm.tsv` | Step 1 |
| `in_kegg_tax` | `kegg_tax.tsv` | Step 1 |
| `in_env` | `metadata_2024_11_11.xlsx` | Physico-chemical data per station |
| `in_ko_path`, `in_path_names`, `in_ko_info` | `ko_path_tbl.rds`, `pathway_names_tbl.rds`, `ko_info_tbl.rds` | KEGG pathway membership and names, retrieved with KEGGREST on [date]. If the first two are missing and KEGGREST is installed, the script downloads them again |
| `in_ncyc`, `in_ncyc_sheet` | `Resultado_FC_only_PN_plus_FV.xlsx` | Abundance of NCycDB gene families per station |
| `in_ncyc_sets` | `Ncyccompleteness.xlsx` | Gene families of each nitrogen-cycle process |

### Structure

| Part | Analysis | Manuscript |
| --- | --- | --- |
| 0 | Configuration, sampling design, palettes, themes, helper functions | – |
| A | Taxonomy from mOTUs: composition, detection limit, Hill diversity at equal coverage, Bray–Curtis ordination, PERMANOVA and PERMDISP | Figure 2, Table 1, Figure S1 |
| B | KEGG orthologues: pathway categories, richness and Hill diversity, sample coverage, Bray–Curtis ordination, PERMANOVA and PERMDISP | Figure 3 |
| C | Comparison of fractions at the four paired stations: rarefaction, shared and exclusive KOs, over-representation analysis, paired log2 ratios per pathway | Figure 4 |
| D | Nitrogen cycle: gene-family richness, nutrients, abundance per process, gene families detected per process | Figure 5 |

### Output

Everything is written to the folder set in `CFG$dir_out` (default `outputs/`).

| File | Content |
| --- | --- |
| `Figure2`, `Figure3`, `Figure4`, `Figure5`, `FigureS1_phylum_composition` (`.jpg` and `.pdf`) | Figures, 600 dpi |
| `Table1_phylum_composition.tsv` | Table 1 |
| `TableS_*.tsv` | Supplementary tables: KEGG categories by fraction, KO diversity and coverage, KO partition by zone, pathway and KO ratios, leading KOs of Figure 4D, over-representation analysis, nitrogen gene families |
| `alpha_taxonomic_hill_q1.tsv` | Hill number of order 1 per sample, with confidence limits |
| `check_*.tsv` | Supporting checks: detection limit, paired dissimilarity, sensitivity to typical contaminant genera, lineage of ORFs |
| `numbers_for_manuscript.tsv` | Every statistic quoted in the text, with its label |
| `sessionInfo.txt` | R and package versions of the run |

## Reproducibility notes

- Random seeds are fixed (`CFG$seed`, default 123) immediately before each permutation test and rarefaction.
- PERMANOVA uses 999 permutations for mOTUs and 9,999 for KOs; PERMDISP uses 999.
- KO presence at the paired stations is assessed after a single rarefaction of the eight samples to the depth of the shallowest one.
- KEGG changes over time. The `.rds` pathway tables are the frozen version used in the thesis.
