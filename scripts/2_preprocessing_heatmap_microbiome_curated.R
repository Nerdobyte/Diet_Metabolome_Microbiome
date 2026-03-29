rm(list=ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(ComplexHeatmap)
  library(circlize)
  library(broom)
})

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
    "to scripts/, or run via Rscript on this file.\ngetwd(): ",
    getwd()
  )
}
repo_root <- resolve_project_root()
input_dir <- file.path(repo_root, "input")
output_dir <- file.path(repo_root, "output")
fig_dir <- file.path(output_dir, "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

require_file <- function(path) {
  if (!file.exists(path)) stop("Required file not found: ", path)
  path
}

metabo_matrix <- readr::read_csv(
  require_file(file.path(output_dir, "output_for_missforest_T22pt1_ALLmol_nodup.csv")),
  na = c("", "NA")
)
if (ncol(metabo_matrix) > 0 && colnames(metabo_matrix)[1] %in% c("X", "...1", "")) {
  metabo_matrix <- metabo_matrix[, -1]
}
metabo_matrix <- janitor::clean_names(metabo_matrix)

df_long <- metabo_matrix %>%
  pivot_longer(
    cols = -c(patient_id, diet, location, timepoint),
    names_to = "molecule",
    values_to = "value"
  ) %>%
  mutate(
    patient_id = as.factor(patient_id),
    diet = as.factor(diet),
    location = as.factor(location),
    timepoint = as.numeric(timepoint)
  )

df_long_filtered <- df_long %>% filter(!(patient_id %in% c(11, 12)))

pid_map <- c(`3009`=1, `9793`=2, `1918`=3, `8520`=4, `2453`=5, `5675`=6, `2588`=7, `7901`=8)

df_ileal <- read.csv(require_file(file.path(input_dir, "All_Ileal_PYY_GLP1.csv")), sep = ",")
colnames(df_ileal) <- gsub("X", "", colnames(df_ileal))
pid_chr <- as.character(df_ileal$PID)
mapped_pid <- unname(pid_map[pid_chr])
df_ileal$PID <- ifelse(!is.na(mapped_pid), mapped_pid, as.numeric(pid_chr))
df_ileal$Location <- "Ileal"

df_gasduo <- read.csv(require_file(file.path(input_dir, "All_GasDuo_PYY_GLP1_GIP.csv")), sep = ",")
colnames(df_gasduo) <- gsub("X", "", colnames(df_gasduo))

df_gastric <- df_gasduo; df_gastric$Location <- "Gastric"
df_duodenal <- df_gasduo; df_duodenal$Location <- "Duodenal"
df_gasduo <- dplyr::bind_rows(df_gastric, df_duodenal)

df_gasduo_long <- df_gasduo %>%
  pivot_longer(cols = c(`0`,`15`,`30`,`45`,`60`,`90`,`120`,`150`,`180`), names_to = "Timepoint", values_to = "Concentration") %>%
  mutate(Timepoint = as.numeric(Timepoint))

df_ileal_long <- df_ileal %>%
  pivot_longer(cols = c(`0`,`60`,`120`,`180`,`240`,`300`,`360`,`420`,`480`), names_to = "Timepoint", values_to = "Concentration") %>%
  mutate(Timepoint = as.numeric(Timepoint))

df_hormones <- bind_rows(df_gasduo_long, df_ileal_long)

df_long_filtered$unit <- "mmol/L"
df_hormones$unit <- "pmol/L"
df_hormones <- janitor::clean_names(df_hormones)
colnames(df_hormones)[3] <- "patient_id"
colnames(df_hormones)[6] <- "value"
df_hormones <- df_hormones %>%
  mutate(
    patient_id = as.factor(patient_id),
    diet = as.factor(recode(as.factor(diet), "BC" = "0", "SC" = "1", "CC" = "2"))
  )

df_combined <- bind_rows(df_long_filtered, df_hormones)
drop_idx <- intersect(c(7, 8), seq_len(ncol(df_combined)))
if (length(drop_idx) > 0) {
  df_combined <- df_combined[, -drop_idx, drop = FALSE]
}
df_combined <- df_combined[df_combined$location != "Oral", ]
df_combined <- na.omit(df_combined)
df_combined$molecule <- gsub("-", "_", df_combined$molecule)

write.csv(df_combined, file.path(output_dir, "df_combined_nonsem.csv"), row.names = FALSE)

df_wide <- df_combined %>%
  pivot_wider(names_from = molecule, values_from = value) %>%
  arrange(patient_id, location, timepoint)

if (ncol(df_wide) >= 32) {
  colnames(df_wide)[16] <- "stachyose + raffinose"
  colnames(df_wide)[22] <- "alpha-glucose + alpha-maltose"
  colnames(df_wide)[23] <- "beta_glucose"
  colnames(df_wide)[24] <- "beta_maltose"
  colnames(df_wide)[31] <- "ciceritol"
  colnames(df_wide)[32] <- "sucrose"
}

molecule_names <- colnames(df_wide)[5:35]
df_wide <- df_wide %>% mutate(diet = recode(diet, "0" = "BC", "1" = "SC", "2" = "CC"))
df_wide <- df_wide[rowSums(is.na(df_wide)) <= 1, ]

df_long_heat <- df_wide %>%
  mutate(diet_time = paste(diet, timepoint, sep = "_T")) %>%
  select(location, diet_time, all_of(molecule_names)) %>%
  pivot_longer(cols = all_of(molecule_names), names_to = "molecule", values_to = "value") %>%
  group_by(molecule) %>%
  mutate(scaled_value = scale(value)[, 1]) %>%
  ungroup()

make_heatmap <- function(loc) {
  df_sub <- df_long_heat %>%
    filter(location == loc) %>%
    select(diet_time, molecule, scaled_value) %>%
    pivot_wider(names_from = diet_time, values_from = scaled_value, values_fn = mean)

  mat <- df_sub %>% tibble::column_to_rownames("molecule") %>% as.matrix()
  mat <- mat[rowSums(is.na(mat)) < ncol(mat), ]

  col_names <- colnames(mat)
  diet <- sub("_T.*", "", col_names)
  timepoint <- as.integer(sub(".*_T", "", col_names))

  ha_col <- HeatmapAnnotation(
    Timepoint = timepoint,
    Diet = diet,
    col = list(
      Timepoint = colorRamp2(c(0, 480), c("white", "black")),
      Diet = c("BC" = "#e75480", "SC" = "black", "CC" = "#6db4b4")
    ),
    annotation_name_side = "right"
  )

  Heatmap(
    mat,
    name = paste0("Z-score_", loc),
    top_annotation = ha_col,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_column_names = FALSE,
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = 10),
    column_split = factor(diet, levels = c("BC", "SC", "CC")),
    col = colorRamp2(c(-2, 0, 2), c("blue", "yellow", "red")),
    na_col = "grey30",
    column_title = loc,
    column_title_gp = gpar(fontsize = 12, fontface = "bold")
  )
}

ht1 <- make_heatmap("Gastric")
ht2 <- make_heatmap("Duodenal")
ht3 <- make_heatmap("Ileal")
pdf(file.path(fig_dir, "location_heatmaps.pdf"), width = 12, height = 14)
draw(ht1, heatmap_legend_side = "right", annotation_legend_side = "right")
draw(ht2, heatmap_legend_side = "right", annotation_legend_side = "right")
draw(ht3, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

df_wide_ileal <- df_combined %>%
  pivot_wider(names_from = molecule, values_from = value) %>%
  arrange(patient_id, location, timepoint) %>%
  subset(location %in% "Ileal")

if ("GIP" %in% colnames(df_wide_ileal)) {
  df_wide_ileal <- df_wide_ileal[, colnames(df_wide_ileal) != "GIP", drop = FALSE]
}
df_wide_ileal <- na.omit(df_wide_ileal)

num_cols <- which(vapply(df_wide_ileal, is.numeric, logical(1)))
num_cols <- num_cols[num_cols > 4]
if (length(num_cols) > 0) df_wide_ileal[, num_cols] <- log1p(df_wide_ileal[, num_cols])

micro <- read.csv(require_file(file.path(input_dir, "Taxa_MGS.matL7b.txt")), sep = "	")
micro[, 2:ncol(micro)] <- log1p(micro[, 2:ncol(micro)])
micro <- micro[-1, ]
colnames(micro)[1] <- "taxon"

micro_long <- micro %>%
  pivot_longer(cols = -taxon, names_to = "sample_id", values_to = "abundance")

linkage <- read.csv(require_file(file.path(input_dir, "Link_microbiome_metabolome_samples.csv")), sep = ",")
linkage <- linkage[-c(138:190), c(1:3)] %>%
  separate(Locator_ICL, into = c("patient", "diet", "timepoint"), sep = "_", convert = TRUE)
colnames(linkage)[4] <- "sample_id"

micro_long <- left_join(micro_long, linkage, by = "sample_id") %>% na.omit()
colnames(micro_long)[4] <- "patient_id"
micro_long$patient_id <- factor(micro_long$patient_id, levels = c(1:8))
micro_long$diet <- factor(micro_long$diet)

df_merged <- micro_long %>% left_join(df_wide_ileal, by = c("patient_id", "diet", "timepoint"))

cor_glp1 <- df_merged %>%
  group_by(timepoint) %>%
  group_map(~ broom::tidy(cor.test(.x$abundance, .x$GLP_1, method = "spearman")), .keep = TRUE) %>%
  bind_rows()

lm_glp1 <- df_merged %>%
  group_by(timepoint, diet) %>%
  group_map(~ broom::tidy(lm(GLP_1 ~ abundance, data = .x)) %>%
              mutate(timepoint = unique(.x$timepoint), diet = unique(.x$diet)), .keep = TRUE) %>%
  bind_rows() %>%
  filter(term == "abundance") %>%
  select(timepoint, diet, estimate, p.value)

write.csv(df_merged, file.path(output_dir, "df_merged_nonsem.csv"), row.names = FALSE)
write.csv(cor_glp1, file.path(output_dir, "cor_glp1_by_timepoint.csv"), row.names = FALSE)
write.csv(lm_glp1, file.path(output_dir, "lm_glp1_by_timepoint_diet.csv"), row.names = FALSE)
