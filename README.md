# Diet_Metabolome_Microbiome Pipeline

Curated R pipeline for gut metabolite imputation (MissForest), wide exports (T21/T22), hormone and microbiome integration, exploratory heatmaps, SEM-style figures, and box plots with peak summaries.

## Layout (project package)

The runnable workflow expects this **three-folder layout** next to each other:

| Folder | Role |
|--------|------|
| `scripts/` | Numbered `*_curated.R` scripts (run in order). |
| `input/` | CSV/TSV inputs you supply (see below). |
| `output/` | All generated tables, figures, and subfolders (created if missing). |

In this repository the package lives at:

`Scripts_Curated/scripts/`, `Scripts_Curated/input/`, `Scripts_Curated/output/`

Scripts **do not hardcode** the parent folder name. They find the project root by locating a directory that contains **`input/`** (and assume `scripts/` and `output/` live beside it). Typical use:

- Run with `Rscript path/to/scripts/1_....R`, or  
- Set the R working directory to the project root (folder that contains `input/`) or to `scripts/`.

## R packages

Install as needed (versions will affect MissForest numerics; pin for full reproducibility):

- **Script 1:** `missForest`, `tidyverse`, `ggfortify`
- **Script 2:** `tidyverse`, `janitor`, `ComplexHeatmap`, `circlize`, `broom`
- **Script 3:** `tidyverse`, `janitor`, `lavaan`, `patchwork`
- **Script 4:** `tidyverse`, `ggpubr`, `janitor`

## Required inputs (`input/`)

| File | Used by |
|------|---------|
| `input_for_missforest_T20_ALLmol_nodup_20250816.csv` | 1 (primary metabolite matrix) |
| `metadata_25072024_for_use_with_metabolites_matrix.csv` | 1 (Halla merge) |
| `All_Ileal_PYY_GLP1.csv` | 2 |
| `All_GasDuo_PYY_GLP1_GIP.csv` | 2 |
| `Taxa_MGS.matL7b.txt` | 2 |
| `Link_microbiome_metabolome_samples.csv` | 2 |

Script **2-4** expect script **1** to have written `output/output_for_missforest_T22pt1_ALLmol_nodup.csv` (and script 1 also produces other T21/T22 files used downstream).

## Run order

From the repository root (adjust path if your clone differs):

```bash
Rscript Scripts_Curated/scripts/1_missforest_imputation_and_exports_curated.R
Rscript Scripts_Curated/scripts/2_preprocessing_heatmap_microbiome_curated.R
Rscript Scripts_Curated/scripts/3_sem_taxon_hormone_heatmap_curated.R
Rscript Scripts_Curated/scripts/4_boxwhisker_plots_and_peak_stats_curated.R
```

### Script 1: imputation and exports

- Reads T20, drops metadata columns, keeps **non-oral** rows for MissForest; **oral** rows are appended again for T22-wide outputs.
- **Per-metabolite** MissForest: each fit uses only `Diet`, `PatientID`, `Location`, `Timepoint`, and **one** metabolite column (no mixing of analytes in a single forest).
- **Bile acids** (TCDCA, GCDCA, TDCA, GDCA, TCA, GCA): trained on **Duodenal + Ileal** only; gastric cells stay as in the input.
- Optional: `MISSFOREST_SEED` (default `20250816` if unset in the environment).

**Main writes under `output/`:**

- `output_for_missforest_T21pt1_ALLmol_nodup.csv` — imputed non-oral wide  
- `output_for_missforest_T22pt1_ALLmol_nodup.csv` — T21 + oral rows  
- `output_for_missforest_T21pt2_ALLmol_nodup.csv`, `output_for_missforest_T22pt2_ALLmol_nodup.csv` — patient-spread (pt2) forms  
- `figures/pca_*.png` — PCA on T21 metabolite block (after median fill for PCA only)  
- `aa_inputs_frameshiftcorrected/aaNN_*.csv` — per-metabolite tables for Origin-style workflows  
- `Halla_20240623/GUT_halla_D*.csv` — diet-stratified matrices for Halla-style use  

### Script 2: preprocessing, heatmaps, microbiome merge

- Combines T22 metabolites (excludes patients 11–12, drops Oral location in the combined long table) with hormone data; writes long combined table.
- Builds **ComplexHeatmap** PDF: `figures/location_heatmaps.pdf` (Gastric, Duodenal, Ileal pages; z-scored values by molecule).
- Merges ileal-wide data with taxon abundances; exports correlation and simple LM summaries for GLP-1 vs abundance.

**Main writes:**

- `df_combined_nonsem.csv` — long metabolites + hormones  
- `df_merged_nonsem.csv` — microbiome × ileal metabolite/hormone merge  
- `cor_glp1_by_timepoint.csv` — Spearman GLP-1 vs abundance by timepoint  
- `lm_glp1_by_timepoint_diet.csv` — `GLP_1 ~ abundance` slope by timepoint and diet  
- `figures/location_heatmaps.pdf`  

### Script 3: SEM / hormone panels

- Reads `output_for_missforest_T22pt1_ALLmol_nodup.csv` and runs lavaan-driven summaries; saves hormone-related plots and a combined taxon–hormone–metabolite heatmap figure under `output/figures/` (e.g. `hormone_*.png`, `hormone_panel_grid.png`, `p_combined_taxon_hormone_metabolite_heatmap.png`).

### Script 4: box plots and peak timepoints

- Reads `output_for_missforest_T22pt1_ALLmol_nodup.csv`, pivots long, and saves PNGs under `output/figures/` in subfolders such as `boxwhisker_basic/`, `boxwhisker_clean_all/`, and `peak_timepoints/` (one set per metabolite column).

## Reproducibility note

MissForest output depends on **R version**, **missForest** and **ranger/randomForest** versions, and **RNG**. The same code and seed can still differ across machines if packages differ. For strict replication, record `sessionInfo()` and package versions after a successful run.

## Publication
"Dietary Modulation of The Human Small Intestine Metabolome and Microbiome"
_(**Nadia Fernandes, Mingzhu Cai, Frederick J Warren,** Anthony Duncan, Hannah Harris, Jose Ivan Serrano Contreras, Katarzyna Sidorczuk, Andres Bernal, Andres Castillo, Katerina Petropoulou, Dominic Blunt, Natalia Perez-Moral, Isabel Garcia-Perez, Cathrina Edwards, Falk Hildebrand, Elaine Holmes, Julien Wist, Gary Frost)_

