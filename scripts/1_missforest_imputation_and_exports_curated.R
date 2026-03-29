rm(list = ls())

suppressPackageStartupMessages({
  library(missForest)
  library(tidyverse)
  library(ggfortify)
})

# Project layout: {project_root}/scripts/*.R, {project_root}/input/, {project_root}/output/
resolve_project_root <- function() {
  has_input_dir <- function(d) {
    d <- normalizePath(d, winslash = "/", mustWork = FALSE)
    dir.exists(file.path(d, "input"))
  }
  aa <- commandArgs(trailingOnly = FALSE)
  fa <- grep("^--file=", aa, value = TRUE)
  if (length(fa)) {
    script_dir <- dirname(normalizePath(sub("^--file=", "", fa[1]), winslash = "/"))
    if (basename(script_dir) == "scripts") {
      candidate <- dirname(script_dir)
      if (has_input_dir(candidate)) {
        return(candidate)
      }
    }
    if (has_input_dir(script_dir)) {
      return(script_dir)
    }
    stop(
      "Could not find project root (folder containing input/ next to scripts/). Script directory: ",
      script_dir
    )
  }
  wd <- normalizePath(getwd(), winslash = "/")
  if (basename(wd) == "scripts") {
    candidate <- dirname(wd)
    if (has_input_dir(candidate)) {
      return(candidate)
    }
  }
  if (has_input_dir(wd)) {
    return(wd)
  }
  subs <- list.dirs(wd, full.names = TRUE, recursive = FALSE)
  for (candidate in subs) {
    if (has_input_dir(candidate)) {
      return(normalizePath(candidate, winslash = "/"))
    }
  }
  stop(
    "Set working directory to the project folder (contains input/ and scripts/), ",
    "to scripts/, or run: Rscript path/to/scripts/this_script.R\ngetwd(): ",
    getwd()
  )
}
repo_root <- resolve_project_root()
input_dir <- file.path(repo_root, "input")
output_dir <- file.path(repo_root, "output")
fig_dir <- file.path(output_dir, "figures")

if (!dir.exists(input_dir)) stop("Input folder does not exist: ", input_dir)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "aa_inputs_frameshiftcorrected"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "Halla_20240623"), recursive = TRUE, showWarnings = FALSE)

safe_read_csv <- function(file_name, ...) {
  in_file <- file.path(input_dir, basename(file_name))
  if (!file.exists(in_file)) stop("Input file not found: ", in_file)
  read.csv(in_file, ...)
}

safe_write_csv <- function(data, file_name, ...) {
  out_file <- file.path(output_dir, basename(file_name))
  write.csv(data, out_file, ...)
  message("Wrote: ", out_file)
}

# Imputation: one missForest call per metabolite. Each run uses exactly these covariates as predictors:
#   Diet, PatientID, Location, Timepoint — plus the single metabolite column being imputed (no other molecules).
# Diet / PatientID / Location are factors; Timepoint is numeric (minutes).
# Six bile acids have no gastric signal in this study: imputation is trained on Duodenal + Ileal rows only; Gastric cells unchanged.
# Env MISSFOREST_SEED (default 20250816) for reproducibility; parallelize = "no" for stable RNG.

gut_raw <- safe_read_csv("input_for_missforest_T20_ALLmol_nodup_20250816.csv", header = TRUE, sep = ",")
gut_raw <- gut_raw[, -c(1, 6:19)]

# Covariates passed into missForest for every molecule (explicit; all four are always included).
mf_covariates <- c("Diet", "PatientID", "Location", "Timepoint")
# Column order in written wide tables (downstream scripts expect these names).
id_cols <- c("PatientID", "Diet", "Location", "Timepoint")
if (!all(mf_covariates %in% colnames(gut_raw)) || !all(id_cols %in% colnames(gut_raw))) {
  stop(
    "T20 is missing required ID columns. Need: ",
    paste(unique(c(mf_covariates, id_cols)), collapse = ", ")
  )
}
if (!setequal(mf_covariates, id_cols)) {
  stop("Internal error: mf_covariates and id_cols must name the same four columns.")
}

gut_raw$.rowid <- seq_len(nrow(gut_raw))

gut_nonoral <- gut_raw %>% filter(.data$Location != "Oral")
gut_oral <- gut_raw %>% filter(.data$Location == "Oral")

gut_nonoral$Location <- as.factor(gut_nonoral$Location)
gut_nonoral$Diet <- as.factor(gut_nonoral$Diet)
gut_nonoral$PatientID <- as.factor(gut_nonoral$PatientID)
gut_nonoral$Timepoint <- suppressWarnings(as.numeric(as.character(gut_nonoral$Timepoint)))
gut_oral$Timepoint <- suppressWarnings(as.numeric(as.character(gut_oral$Timepoint)))

# Keep all metabolite columns from the non-oral subset; oral rows are added back only for T22 outputs.
metabolite_cols <- colnames(gut_nonoral)[!colnames(gut_nonoral) %in% c(".rowid", id_cols)]
if (length(metabolite_cols) == 0L) stop("No metabolite columns left after exclusions.")

bile_acid_cols <- intersect(
  metabolite_cols,
  c("TCDCA", "GCDCA", "TDCA", "GDCA", "TCA", "GCA")
)

imputation_seed <- suppressWarnings(as.integer(Sys.getenv("MISSFOREST_SEED", unset = "20250816")))
if (is.na(imputation_seed)) imputation_seed <- 20250816L
set.seed(imputation_seed)

# --- Per-metabolite imputation: never pass more than one response column to missForest ---
impute_one_metabolite <- function(mcol) {
  full_vec <- gut_nonoral[[mcol]]
  stopifnot(length(full_vec) == nrow(gut_nonoral))

  if (all(is.na(full_vec))) {
    return(full_vec)
  }

  if (mcol %in% bile_acid_cols) {
    out <- full_vec
    idx_train <- gut_nonoral$Location %in% c("Duodenal", "Ileal")
    mat_train <- gut_nonoral[idx_train, c(mf_covariates, mcol), drop = FALSE]
    if (all(is.na(mat_train[[mcol]]))) {
      return(out)
    }
    mf <- missForest(mat_train, verbose = FALSE, parallelize = "no")
    if (!(mcol %in% colnames(mf$ximp))) {
      return(out)
    }
    imp <- mf$ximp[[mcol]]
    stopifnot(length(imp) == sum(idx_train))
    out[idx_train] <- imp
    return(out)
  }

  mat <- gut_nonoral[, c(mf_covariates, mcol), drop = FALSE]
  mf <- missForest(mat, verbose = FALSE, parallelize = "no")
  if (!(mcol %in% colnames(mf$ximp))) {
    return(full_vec)
  }
  imp <- mf$ximp[[mcol]]
  stopifnot(length(imp) == nrow(gut_nonoral))
  imp
}

imputed_list <- vector("list", length(metabolite_cols))
names(imputed_list) <- metabolite_cols
for (mcol in metabolite_cols) {
  imputed_list[[mcol]] <- impute_one_metabolite(mcol)
}

# Assemble T21pt1: preserve row order and column names exactly.
gut_t21pt1 <- gut_nonoral[, c(".rowid", id_cols), drop = FALSE]
for (mcol in metabolite_cols) {
  gut_t21pt1[[mcol]] <- imputed_list[[mcol]]
}

safe_write_csv(gut_t21pt1 %>% select(-.rowid), "output_for_missforest_T21pt1_ALLmol_nodup.csv", row.names = FALSE)

gut_oral_ordered <- gut_oral[, c(".rowid", id_cols, metabolite_cols), drop = FALSE]
gut_t21pt1 <- gut_t21pt1 %>% mutate(across(all_of(id_cols), as.character))
gut_oral_ordered <- gut_oral_ordered %>% mutate(across(all_of(id_cols), as.character))
gut_t22pt1 <- bind_rows(gut_t21pt1, gut_oral_ordered) %>%
  arrange(.rowid) %>%
  select(.rowid, all_of(id_cols), all_of(metabolite_cols))
safe_write_csv(gut_t22pt1 %>% select(-.rowid), "output_for_missforest_T22pt1_ALLmol_nodup.csv", row.names = FALSE)


gut_t21 <- gut_t21pt1 %>% select(-.rowid)
gut_t22 <- gut_t22pt1 %>% select(-.rowid)

to_spread <- function(df_in) {
  df_in %>%
    mutate(Location_Timepoint = str_c(.data$Location, "_", .data$Timepoint)) %>%
    relocate(Location_Timepoint, .after = Timepoint) %>%
    select(-Location) %>%
    pivot_longer(cols = all_of(metabolite_cols), names_to = "Molecule", values_to = "Conc") %>%
    pivot_wider(names_from = PatientID, values_from = Conc)
}

gut_spread_t21 <- to_spread(gut_t21)
gut_spread_t22 <- to_spread(gut_t22)

safe_write_csv(gut_spread_t21, "output_for_missforest_T21pt2_ALLmol_nodup.csv", row.names = FALSE)
safe_write_csv(gut_spread_t22, "output_for_missforest_T22pt2_ALLmol_nodup.csv", row.names = FALSE)

gut_imputed <- gut_t22
gut_spread <- gut_spread_t22
gut_imputed$Location_Timepoint <- str_c(gut_imputed$Location, "_", gut_imputed$Timepoint)
gut_long <- gut_imputed %>%
  relocate(Location_Timepoint, .after = Timepoint) %>%
  select(-Location) %>%
  pivot_longer(cols = all_of(metabolite_cols), names_to = "Molecule", values_to = "Conc")


gut_imputed_old <- gut_imputed
gut_imputed <- read.csv(file.path(output_dir, "output_for_missforest_T21pt1_ALLmol_nodup.csv"), header = TRUE, sep = ",")
# PCA on all imputed metabolite columns from T21pt1; coerce to numeric and median-fill NA/Inf for prcomp only.
df1 <- gut_imputed[, metabolite_cols, drop = FALSE]
df1 <- as.data.frame(lapply(df1, function(v) {
  v <- suppressWarnings(as.numeric(v))
  v[!is.finite(v)] <- NA_real_
  med <- stats::median(v, na.rm = TRUE)
  if (is.na(med)) med <- 0
  v[is.na(v)] <- med
  v
}))
sd_ok <- vapply(df1, stats::sd, numeric(1), na.rm = TRUE) > 0
if (!all(sd_ok)) {
  message("PCA: dropping ", sum(!sd_ok), " zero-variance metabolite column(s).")
}
df1 <- df1[, sd_ok, drop = FALSE]
pca_res <- prcomp(df1, center = TRUE, scale. = TRUE)

gut_imputed$Diet <- as.factor(gut_imputed$Diet)
gut_imputed$Location <- as.factor(gut_imputed$Location)
gut_imputed$Timepoint <- as.factor(gut_imputed$Timepoint)
gut_imputed$PatientID <- as.factor(gut_imputed$PatientID)
p_patient <- autoplot(pca_res, data = gut_imputed, colour = "PatientID")
p_location <- autoplot(pca_res, data = gut_imputed, colour = "Location")
gut_imputed$Diet <- gsub("0", "BC", gut_imputed$Diet)
gut_imputed$Diet <- gsub("1", "SC", gut_imputed$Diet)
gut_imputed$Diet <- gsub("2", "CC", gut_imputed$Diet)
p_diet <- autoplot(pca_res, data = gut_imputed, colour = "Diet")
p_timepoint <- autoplot(pca_res, data = gut_imputed, colour = "Timepoint")

ggsave(file.path(fig_dir, "pca_patientid.png"), p_patient, width = 7, height = 5, dpi = 300)
ggsave(file.path(fig_dir, "pca_location.png"), p_location, width = 7, height = 5, dpi = 300)
ggsave(file.path(fig_dir, "pca_diet.png"), p_diet, width = 7, height = 5, dpi = 300)
ggsave(file.path(fig_dir, "pca_timepoint.png"), p_timepoint, width = 7, height = 5, dpi = 300)

gut_imputed <- gut_imputed_old


to_originlab_matrix <- function(molecule_name) {
  aa_data <- gut_spread %>%
    filter(Molecule == molecule_name) %>%
    mutate(means = rowMeans(across(where(is.numeric)), na.rm = TRUE)) %>%
    select(Location_Timepoint, Diet, means) %>%
    separate(Location_Timepoint, into = c("Location", "Timepoint"), sep = "_(?=[^_]+$)") %>%
    mutate(Diet_Location = str_c(Diet, "_", Location)) %>%
    select(Timepoint, Diet_Location, means) %>%
    pivot_wider(names_from = Diet_Location, values_from = means)

  aa_data
}

for (idx in seq_along(metabolite_cols)) {
  molecule <- metabolite_cols[idx]
  aa_data <- to_originlab_matrix(molecule)
  file_name <- sprintf("aa%02d_data_frameshiftcorrected.csv", idx)
  out_path <- file.path(output_dir, "aa_inputs_frameshiftcorrected", file_name)
  write.csv(aa_data, out_path, row.names = FALSE, na = "")
}


gut_halla <- gut_spread %>%
  mutate(Location = str_extract(Location_Timepoint, "^[^_]+")) %>%
  filter(Location %in% c("Ileal", "Oral")) %>%
  separate(Location_Timepoint, into = c("Location", "Timepoint"), sep = "_(?=[^_]+$)") %>%
  pivot_longer(cols = where(is.numeric), names_to = "PatientID", values_to = "Conc") %>%
  pivot_wider(names_from = Molecule, values_from = Conc) %>%
  mutate(Locator_ICL = paste0(PatientID, "_", Diet, "_", Timepoint))

meta <- safe_read_csv("metadata_25072024_for_use_with_metabolites_matrix.csv", sep = ",", header = TRUE)
gut_halla <- merge(gut_halla, meta, by = "Locator_ICL", all = TRUE)
gut_halla <- gut_halla[!is.na(gut_halla$sample), ]
rownames(gut_halla) <- gut_halla$sample

for (diet_code in c("0", "1", "2")) {
  df <- gut_halla[grepl(diet_code, gut_halla$Diet), ]
  df <- df[, setdiff(colnames(df), "Diet"), drop = FALSE]
  out <- t(df)
  out_name <- sprintf("GUT_halla_D%s.csv", diet_code)
  out_path <- file.path(output_dir, "Halla_20240623", out_name)
  write.csv(out, out_path)
}

message("Script completed successfully.")
