# Diet-Metabolomics-Microbiome Pipeline

Curated R pipeline for gut metabolite imputation (MissForest), wide exports (T21/T22), hormone and microbiome integration, exploratory heatmaps, SEM-style figures, and box plots with peak summaries.

## Layout (project package)

The runnable workflow expects this **three-folder layout** next to each other:

| Folder | Role |
|--------|------|
| `scripts/` | Numbered `*_curated.R` scripts (run in order). |
| `input/` | CSV/TSV inputs you supply (see below). |
| `output/` | All generated tables, figures, and subfolders (created if missing). |

### Path resolution

- Scripts **do not** embed a fixed parent folder name (e.g. no hardcoded machine path).
- The **project root** is resolved as the directory that contains **`input/`**:
  - `Rscript …/scripts/<script>.R` → parent of `scripts/` if that parent has `input/`.
  - Or `getwd()` is the project root (has `input/`), or `getwd()` is `scripts/` (then parent is used), or a direct child of `getwd()` contains `input/`.

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

Scripts **2–4** expect script **1** to have written `output/output_for_missforest_T22pt1_ALLmol_nodup.csv`.

## Run order

From the repository root (adjust if your clone path differs):

```bash
Rscript scripts/1_missforest_imputation_and_exports_curated.R
Rscript scripts/2_preprocessing_heatmap_microbiome_curated.R
Rscript scripts/3_sem_taxon_hormone_heatmap_curated.R
Rscript scripts/4_boxwhisker_plots_and_peak_stats_curated.R
```

### Script 1: imputation and exports

- Reads T20, drops index + metadata columns (`1` and `6:19` from the wide table).
- **Non-oral** rows only are passed to MissForest; **oral** rows from T20 are bound back for **T22** wide outputs (not imputed in script 1).
- **All** metabolite columns present in the non-oral table are imputed (no separate drop list of “oral-only” analytes).
- **Per-metabolite** MissForest: each fit uses **`Diet`, `PatientID`, `Location`, `Timepoint`** (predictors) plus **one** metabolite column only — no mixing of molecules in a single forest.
- **Wide CSV column order:** `PatientID`, `Diet`, `Location`, `Timepoint`, then metabolites (same names as T20 after the column drop).
- **Bile acids** (TCDCA, GCDCA, TDCA, GDCA, TCA, GCA): trained on **Duodenal + Ileal** only; **Gastric** cells stay as in T20.
- Environment: **`MISSFOREST_SEED`** (default `20250816` if unset). Uses `missForest(..., parallelize = "no")`.

**Main writes under `output/`:**

- `output_for_missforest_T21pt1_ALLmol_nodup.csv` — imputed non-oral wide  
- `output_for_missforest_T22pt1_ALLmol_nodup.csv` — T21 + oral rows  
- `output_for_missforest_T21pt2_ALLmol_nodup.csv`, `output_for_missforest_T22pt2_ALLmol_nodup.csv` — patient-spread (pt2) forms  
- `figures/pca_*.png` — PCA on T21 metabolite columns after numeric coercion + median imputation for missing/non-finite values; **columns with zero variance after that step are dropped** before `prcomp` (console message when that happens)  
- `aa_inputs_frameshiftcorrected/aaNN_*.csv` — per-metabolite tables for Origin-style workflows  
- `Halla_20240623/GUT_halla_D*.csv` — diet-stratified matrices for Halla-style use  

### Script 2: preprocessing, heatmaps, microbiome merge

- Combines T22 metabolites (excludes patients 11–12, drops **Oral** in the combined long table) with hormone data.
- **ComplexHeatmap** PDF: `figures/location_heatmaps.pdf` (Gastric, Duodenal, Ileal; z-scored by molecule).
- Merges ileal-wide data with taxon abundances; exports GLP-1 vs abundance correlation and simple LM slices.

**Main writes:**

- `df_combined_nonsem.csv` — long metabolites + hormones  
- `df_merged_nonsem.csv` — microbiome × ileal metabolite/hormone merge  
- `cor_glp1_by_timepoint.csv` — Spearman GLP-1 vs abundance by timepoint  
- `lm_glp1_by_timepoint_diet.csv` — `GLP_1 ~ abundance` coefficient by timepoint and diet  
- `figures/location_heatmaps.pdf`  

### Script 3: SEM/ hormone panels

- Reads `output/output_for_missforest_T22pt1_ALLmol_nodup.csv`, lavaan-driven summaries.
- Figures under `output/figures/`, e.g. `hormone_*.png`, `hormone_panel_grid.png`, `p_combined_taxon_hormone_metabolite_heatmap.png`.

### Script 4: box plots and peak timepoints

- Reads `output/output_for_missforest_T22pt1_ALLmol_nodup.csv` and builds long-format plotting data (`janitor::clean_names` → `patient_id`, `diet`, `location`, `timepoint`).
- **Writes only two PNG families** (no diet×timepoint faceted “clean all” grid):
  - **`figures/boxwhisker_basic/`** — `box_basic_<metabolite>.png`: timepoint on x, concentration on y, diet fill, facets by location.
  - **`figures/peak_timepoints/`** — `peak_<metabolite>.png`: per patient/location/diet, timepoint of max concentration (non-oral); optional significance brackets via `ggpubr`.
- **`boxwhisker_clean_all/` is not generated** by the current script (removed to shorten run time and output size).

## Reproducibility note

MissForest output depends on **R version**, **missForest** and **ranger/randomForest** versions, and **RNG**. The same code and seed can still differ across machines if packages differ. For strict replication, record `sessionInfo()` and package versions after a successful run.

## Publication

**Dietary Modulation of The Human Small Intestine Metabolome and Microbiome**

*Nadia Fernandes, Mingzhu Cai, Frederick J Warren, Anthony Duncan, Hannah Harris, Jose Ivan Serrano Contreras, Katarzyna Sidorczuk, Andres Bernal, Andres Castillo, Katerina Petropoulou, Dominic Blunt, Natalia Perez-Moral, Isabel Garcia-Perez, Cathrina Edwards, Falk Hildebrand, Elaine Holmes, Julien Wist, Gary Frost*
