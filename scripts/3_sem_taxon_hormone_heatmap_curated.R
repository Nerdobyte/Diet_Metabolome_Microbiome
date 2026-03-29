rm(list=ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(lavaan)
  library(patchwork)
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

metabo_matrix = read_csv(file.path(output_dir, "output_for_missforest_T22pt1_ALLmol_nodup.csv"), na = c("", "NA"))
if (ncol(metabo_matrix) > 0 && colnames(metabo_matrix)[1] %in% c("X", "...1", "")) {
  metabo_matrix <- metabo_matrix[, -1]
}
glimpse(metabo_matrix)
metabo_matrix <- janitor::clean_names(metabo_matrix)

chickpea_rows <- metabo_matrix %>%
  filter(str_detect(patient_id, regex("baseline_chickpeas_puree", ignore_case = TRUE)))

print(chickpea_rows)

df_long <- metabo_matrix %>%
  pivot_longer(
    cols = -c(patient_id, diet, location, timepoint),
    names_to = "molecule",
    values_to = "value"
  )


df_long <- df_long %>%
  mutate(
    patient_id = as.factor(patient_id),
    diet = as.factor(diet),
    location = as.factor(location),
    timepoint = as.numeric(timepoint)
  )

df_long_filtered <- df_long %>%
  filter(!(patient_id %in% c(11, 12)))
unique(df_long_filtered$patient_id)


df_long_filtered %>%
  filter(!is.na(value)) %>%
  group_by(patient_id, location, timepoint) %>%
  summarise(n_obs = n(), .groups = "drop") %>%
  count(patient_id)


df_ileal <- read.csv(file.path(input_dir, "All_Ileal_PYY_GLP1.csv"),
                     sep=",")
colnames(df_ileal) = gsub("X","",colnames(df_ileal))
df_ileal$PID[df_ileal$PID == 3009] <- 1
df_ileal$PID[df_ileal$PID == 9793] <- 2
df_ileal$PID[df_ileal$PID == 1918] <- 3
df_ileal$PID[df_ileal$PID == 8520] <- 4
df_ileal$PID[df_ileal$PID == 2453] <- 5
df_ileal$PID[df_ileal$PID == 5675] <- 6
df_ileal$PID[df_ileal$PID == 2588] <- 7
df_ileal$PID[df_ileal$PID == 7901] <- 8
df_gasduo <- read.csv(file.path(input_dir, "All_GasDuo_PYY_GLP1_GIP.csv"),
                      sep=",")
colnames(df_gasduo) = gsub("X","",colnames(df_gasduo))

df_gastric <- df_gasduo
df_duodenal <- df_gasduo

df_gastric$Location <- "Gastric"
df_duodenal$Location <- "Duodenal"

df_gasduo <- rbind(df_gastric, df_duodenal)

print(df_gasduo)

df_gasduo_long <- df_gasduo %>%
  pivot_longer(cols = c(`0`, `15`, `30`, `45`, `60`, `90`, `120`, `150`, `180`),
               names_to = "Timepoint", values_to = "Concentration") %>%
  mutate(Timepoint = as.numeric(Timepoint))

df_ileal$Location <- "Ileal"

df_ileal_long <- df_ileal %>%
  pivot_longer(cols = c(`0`, `60`, `120`, `180`, `240`, `300`, `360`, `420`, `480`),
               names_to = "Timepoint", values_to = "Concentration") %>%
  mutate(Timepoint = as.numeric(Timepoint))

df_hormones <- bind_rows(df_gasduo_long, df_ileal_long)

df_hormones <- df_hormones %>%
  mutate(Phase = case_when(
    Timepoint <= 30 ~ "Induction",
    Timepoint > 30 & Timepoint <= 120 ~ "Maintenance",
    Timepoint > 120 ~ "Extension"
  ))

target_combos <- list(
  Duodenal_GIP = c("Duodenal", "GIP"),
  Duodenal_GLP1 = c("Duodenal", "GLP-1"),
  Duodenal_PYY = c("Duodenal", "PYY"),
  Ileal_GLP1 = c("Ileal", "GLP-1"),
  Ileal_PYY = c("Ileal", "PYY")
)

phase_shading <- data.frame(
  Phase = c("Induction", "Maintenance", "Extension"),
  xmin = c(0, 30, 120),
  xmax = c(30, 120, 480),
  fill = c("grey", "lightgoldenrod1", "skyblue2")
)

plot_hormone_location <- function(df, location, hormone, title_name) {
  sub_df <- df %>% filter(Location == location, Molecule == hormone)
  ggplot(sub_df, aes(x = Timepoint, y = Concentration, color = Diet)) +
    geom_rect(
      data = phase_shading,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = Phase),
      inherit.aes = FALSE, alpha = 0.3
    ) +
    stat_summary(fun = mean, geom = "line", aes(group = Diet)) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 5) +
    stat_summary(fun = mean, geom = "point") +
    scale_fill_manual(values = setNames(phase_shading$fill, phase_shading$Phase)) +
    labs(
      title = title_name,
      x = "Time (min)",
      y = paste(hormone, "(pmol/L)")
    ) +
    theme_minimal(base_size = 14)
}

plots <- list()
for (name in names(target_combos)) {
  loc_h <- target_combos[[name]]
  plots[[name]] <- plot_hormone_location(df_hormones, loc_h[1], loc_h[2], title_name = paste(loc_h[1], loc_h[2], sep = " - "))
}

for (name in names(plots)) {
  ggsave(file.path(fig_dir, paste0("hormone_", name, ".png")), plot = plots[[name]], width = 8, height = 5, dpi = 300)
}

hormone_grid <- (plots$Duodenal_GIP | plots$Duodenal_GLP1 | plots$Duodenal_PYY) /
  (plots$Ileal_GLP1 | plots$Ileal_PYY)
ggsave(file.path(fig_dir, "hormone_panel_grid.png"), plot = hormone_grid, width = 16, height = 10, dpi = 300)

rm(df_duodenal, df_gasduo, df_gastric, df_ileal, df_gasduo_long, df_ileal_long)


df_long_filtered$unit <- "mmol/L"
df_hormones$unit <- "pmol/L"

df_hormones = janitor::clean_names(df_hormones)
colnames(df_hormones)[3] = "patient_id"
colnames(df_hormones)[6] = "value"
df_hormones$patient_id = as.factor(df_hormones$patient_id)
df_hormones$diet = as.factor(df_hormones$diet)
df_hormones <- df_hormones %>%
  mutate(diet = recode(diet, "BC" = 0, "SC" = 1, "CC" = 2))
df_hormones$diet <- as.factor(df_hormones$diet)

df_combined <- bind_rows(df_long_filtered, df_hormones)
df_combined = df_combined[,-c(7,8)]
df_combined <- df_combined[df_combined$location != "Oral", ]
df_combined <- na.omit(df_combined)
df_combined$molecule <- gsub("-", "_", df_combined$molecule)

df_wide <- df_combined %>%
  pivot_wider(names_from = molecule, values_from = value) %>%
  arrange(patient_id, location, timepoint)
df_wide = subset(df_wide, location %in% "Ileal")

df_wide = df_wide[,-35]
df_wide <- na.omit(df_wide)

df_wide[ , c(5:34)] = log1p(df_wide[ , c(5:34)])

micro = read.csv(file.path(input_dir, "Taxa_MGS.matL7b.txt"), sep = "\t")
micro[ , c(2:189)] = log1p(micro[ , c(2:189)])
micro = micro[-1,]
colnames(micro)[1] = "taxon"

micro_long <- micro %>%
  pivot_longer(
    cols = -taxon,
    names_to = "sample_id",
    values_to = "abundance"
  )

linkage = read.csv(file.path(input_dir, "Link_microbiome_metabolome_samples.csv"), sep=",")
linkage = linkage[-c(138:190),c(1:3)]

linkage = linkage %>%
  separate(Locator_ICL, into = c("patient", "diet", "timepoint"), sep = "_", convert = TRUE)

colnames(linkage)[4] = "sample_id"

micro_long <- left_join(micro_long, linkage, by = "sample_id")
micro_long <- na.omit(micro_long)
colnames(micro_long)[4] = "patient_id"
micro_long$patient_id = factor(micro_long$patient_id, levels = c(1:8))
micro_long$diet = factor(micro_long$diet)

df_merged <- micro_long %>%
  left_join(df_wide, by = c("patient_id", "diet", "timepoint"))

colnames(df_merged)[20] = "stachyose + raffinose"
colnames(df_merged)[26] = "alpha-glucose + alpha-maltose"
colnames(df_merged)[27] = "beta_glucose" 
colnames(df_merged)[28] = "beta_maltose"
colnames(df_merged)[35] = "ciceritol"
colnames(df_merged)[36] = "sucrose"

total_samples <- n_distinct(df_merged$sample_id)

taxon_sample_counts <- df_merged %>%
  filter(abundance > 0) %>%
  group_by(taxon) %>%
  summarise(samples_present = n_distinct(sample_id))

threshold <- 0.30 * total_samples
taxa_to_keep1 <- taxon_sample_counts %>%
  filter(samples_present >= threshold) %>%
  pull(taxon)

df_filtered1 <- df_merged %>%
  filter(taxon %in% taxa_to_keep1)

taxa_retained1 <- n_distinct(df_filtered1$taxon)
taxa_total <- n_distinct(df_merged$taxon)
pct_taxa_lost1 <- 100 * (1 - taxa_retained1 / taxa_total)

taxon_avg_abundance <- df_merged %>%
  group_by(taxon) %>%
  summarise(mean_abundance = mean(abundance))

taxa_to_keep2 <- taxon_avg_abundance %>%
  filter(mean_abundance > 1) %>%
  pull(taxon)

df_filtered2 <- df_merged %>%
  filter(taxon %in% taxa_to_keep2)

taxa_retained2 <- n_distinct(df_filtered2$taxon)
pct_taxa_lost2 <- 100 * (1 - taxa_retained2 / taxa_total)

df = df_filtered2

taxa <- unique(df$taxon)
metabolites <- colnames(df)[9:36]
diets <- unique(df$diet)


results_taxon_metabolite <- data.frame(
  taxon = character(),
  diet = character(),
  timepoint = character(),
  metabolite = character(),
  estimate = numeric(),
  pvalue = numeric(),
  stringsAsFactors = FALSE
)

for (t in taxa) {
  for (d in diets) {
    for (m in metabolites) {
      df_sub <- df %>% filter(taxon == t, diet == d)
      if (nrow(df_sub) >= 10) {
        model <- paste0(m, " ~ abundance")
        fit <- tryCatch(sem(model, data = df_sub, estimator = "MLM"), error = function(e) NULL)
        if (!is.null(fit)) {
          coef <- standardizedSolution(fit)
          path <- coef %>% filter(lhs == m & rhs == "abundance")
          results_taxon_metabolite <- rbind(results_taxon_metabolite, data.frame(
            taxon = t,
            diet = d,
            metabolite = m,
            estimate = path$est,
            pvalue = path$pvalue
          ))
        }
      }
    }
  }
}

results_taxon_metabolite$padj <- p.adjust(results_taxon_metabolite$pvalue, method = "BH")

results_taxon_metabolite <- results_taxon_metabolite %>% arrange(pvalue)

results_taxon_metabolite_filtered_extreme <- results_taxon_metabolite %>%
  filter(padj < 0.05 & (estimate < -0.5 | estimate > 0.5))


results_metabolite_hormone <- data.frame(
  metabolite = character(),
  diet = character(),
  hormone = character(),
  estimate = numeric(),
  pvalue = numeric(),
  stringsAsFactors = FALSE
)

hormones <- c("PYY", "GLP_1")

for (m in metabolites) {
  for (d in diets) {
    for (h in hormones) {
      df_sub <- df %>% filter(diet == d, !is.na(df[[m]]), !is.na(get(h)))
      if (nrow(df_sub) >= 10) {
        model <- paste0(h, " ~ ", m)
        fit <- tryCatch(sem(model, data = df_sub, estimator = "MLM"), error = function(e) NULL)
        if (!is.null(fit)) {
          coef <- standardizedSolution(fit)
          path <- coef %>% filter(lhs == h & rhs == m)
          results_metabolite_hormone <- rbind(results_metabolite_hormone, data.frame(
            metabolite = m,
            diet = d,
            hormone = h,
            estimate = path$est,
            pvalue = path$pvalue
          ))
        }
      }
    }
  }
}

results_metabolite_hormone$padj <- p.adjust(results_metabolite_hormone$pvalue, method = "BH")

results_metabolite_hormone <- results_metabolite_hormone %>% arrange(pvalue)

results_metabolite_hormone_filtered <- results_metabolite_hormone %>%
  filter(padj < 0.05)


results_hormone_as_taxon <- results_metabolite_hormone_filtered %>%
  rename(taxon = hormone) %>%
  select(taxon, diet, metabolite, estimate, pvalue, padj)

combined_results <- bind_rows(
  results_taxon_metabolite_filtered_extreme %>% select(taxon, diet, metabolite, estimate, pvalue, padj),
  results_hormone_as_taxon
)

combined_results$taxon <- factor(combined_results$taxon, levels = c(sort(unique(results_taxon_metabolite_filtered_extreme$taxon)), c("GLP_1", "PYY")))

p_combined <- ggplot(combined_results, aes(x = metabolite, y = taxon, fill = estimate)) +
  geom_tile() +
  facet_wrap(~ diet) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  labs(
    title = "Taxon and Hormone Associations with Metabolites per Diet",
    x = "Metabolite",
    y = "Taxon / Hormone",
    fill = "Association Estimate"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8)
  )

ggsave(file.path(fig_dir, "p_combined_taxon_hormone_metabolite_heatmap.png"), plot = p_combined, width = 13, height = 8, dpi = 300)


