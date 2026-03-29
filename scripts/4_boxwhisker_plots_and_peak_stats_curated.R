rm(list=ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
  library(janitor)
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

metabo_long <- metabo_matrix %>%
  pivot_longer(
    cols = -c(patient_id, diet, location, timepoint),
    names_to = "metabolite",
    values_to = "concentration"
  ) %>%
  mutate(
    patient_id = as.character(patient_id),
    diet_num = suppressWarnings(as.numeric(as.character(diet))),
    location = factor(location, levels = c("Oral", "Gastric", "Duodenal", "Ileal")),
    timepoint = as.numeric(timepoint)
  ) %>%
  mutate(
    diet = factor(diet_num, levels = c(0, 1, 2)),
    diet_label = case_when(
      diet_num == 0 ~ "BC",
      diet_num == 1 ~ "SC",
      diet_num == 2 ~ "CC",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("BC", "SC", "CC"))
  )

plot_metabolite <- function(met_name) {
  metabo_long %>%
    filter(metabolite == met_name) %>%
    ggplot(aes(x = factor(timepoint), y = concentration, fill = diet)) +
    geom_boxplot(
      outlier.shape = NA,
      position = position_dodge(width = 0.8),
      width = 0.6
    ) +
    geom_jitter(
      aes(color = diet),
      position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
      size = 1.5,
      alpha = 0.6
    ) +
    facet_wrap(~ location, ncol = 2, scales = "free_y") +
    labs(
      title = met_name,
      x = "Time (min)",
      y = "Concentration",
      fill = "Diet",
      color = "Diet"
    ) +
    scale_fill_manual(values = c("0" = "#1b9e77", "1" = "#d95f02", "2" = "#7570b3")) +
    scale_color_manual(values = c("0" = "#1b9e77", "1" = "#d95f02", "2" = "#7570b3")) +
    theme_minimal(base_size = 12) +
    theme(
      strip.text = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

plot_peak_timepoints <- function(met_name) {
  df_peak_tp <- metabo_long %>%
    filter(metabolite == met_name, location != "Oral") %>%
    group_by(patient_id, location, diet_label) %>%
    slice_max(order_by = concentration, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      timepoint_label = ifelse(timepoint == -1, "Oral", as.character(timepoint)),
      timepoint_label = factor(
        timepoint_label,
        levels = c("0", "12", "15", "30", "60", "90", "120", "135", "150",
                   "165", "180", "240", "300", "360", "480")
      ),
      location = factor(location, levels = c("Gastric", "Duodenal", "Ileal")),
      diet_label = factor(diet_label, levels = c("BC", "SC", "CC")),
      loc_diet = factor(paste(location, diet_label, sep = "_"),
                        levels = c("Gastric_BC", "Gastric_SC", "Gastric_CC",
                                   "Duodenal_BC", "Duodenal_SC", "Duodenal_CC",
                                   "Ileal_BC", "Ileal_SC", "Ileal_CC"))
    ) %>%
    mutate(timepoint_num = as.numeric(as.character(timepoint_label)))

  diet_levels_peak <- c("BC", "SC", "CC")

  stats <- df_peak_tp %>%
    group_by(location) %>%
    group_modify(~ {
      dsub <- .x
      diets_here <- intersect(diet_levels_peak, unique(as.character(dsub$diet_label)))
      if (length(diets_here) < 2) return(tibble())
      pairs <- combn(diets_here, 2, simplify = FALSE)
      purrr::map_dfr(pairs, function(pair) {
        d1 <- filter(dsub, as.character(diet_label) == pair[1])
        d2 <- filter(dsub, as.character(diet_label) == pair[2])
        if (nrow(d1) >= 3 & nrow(d2) >= 3) {
          test <- wilcox.test(d1$timepoint_num, d2$timepoint_num, paired = FALSE)
          tibble(
            location = unique(.x$location),
            diet1 = pair[1],
            diet2 = pair[2],
            p = test$p.value
          )
        } else {
          tibble()
        }
      })
    }) %>%
    ungroup()

  stats <- stats %>%
    mutate(
      p.signif = case_when(
        p <= 0.001 ~ "***",
        p <= 0.01 ~ "**",
        p <= 0.05 ~ "*",
        TRUE ~ "ns"
      ),
      group1 = paste(location, diet1, sep = "_"),
      group2 = paste(location, diet2, sep = "_")
    )

  y_pos_df <- df_peak_tp %>%
    group_by(loc_diet) %>%
    summarise(y_max = max(timepoint_num, na.rm = TRUE), .groups = "drop")

  stats <- stats %>%
    left_join(y_pos_df, by = c("group1" = "loc_diet")) %>%
    rename(y1 = y_max) %>%
    left_join(y_pos_df, by = c("group2" = "loc_diet")) %>%
    rename(y2 = y_max) %>%
    group_by(location) %>%
    arrange(p) %>%
    mutate(
      base_y = pmax(y1, y2),
      y.position = base_y + (row_number() * 30)
    ) %>%
    ungroup()

  p <- ggplot(df_peak_tp, aes(x = loc_diet, y = timepoint_num, fill = diet_label)) +
    geom_boxplot(outlier.shape = NA, width = 0.7) +
    geom_jitter(color = "black", width = 0.15, size = 1.8, alpha = 0.7) +
    scale_fill_manual(values = c("BC" = "#F8766D", "SC" = "#609DFF", "CC" = "#00BA39")) +
    scale_y_continuous(
      breaks = c(0, 12, 30, 60, 90, 120, 180, 240, 360, 480),
      labels = c("0", "12", "30", "60", "90", "120", "180", "240", "360", "480")
    ) +
    labs(
      title = paste(met_name),
      x = "Location_Diet",
      y = "Timepoint of Peak Concentration (minutes)"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      legend.position = "none",
      plot.title = element_text(size = 12, face = "bold")
    )

  if (nrow(stats) > 0) {
    p <- p + stat_pvalue_manual(
      stats,
      label = "p.signif",
      xmin = "group1",
      xmax = "group2",
      y.position = "y.position",
      tip.length = 0.01,
      bracket.size = 0.6,
      size = 4
    )
  }

  p
}

safe_met_file <- function(m) gsub("[^A-Za-z0-9._-]+", "_", m)

dir_box_basic <- file.path(fig_dir, "boxwhisker_basic")
dir_peak <- file.path(fig_dir, "peak_timepoints")
dir.create(dir_box_basic, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_peak, recursive = TRUE, showWarnings = FALSE)

metabolites_all <- sort(unique(metabo_long$metabolite))
message("Saving plots for ", length(metabolites_all), " metabolites under ", fig_dir)

for (m in metabolites_all) {
  safe <- safe_met_file(m)
  tryCatch(
    ggsave(
      file.path(dir_box_basic, paste0("box_basic_", safe, ".png")),
      plot = plot_metabolite(m),
      width = 12,
      height = 10,
      dpi = 300
    ),
    error = function(e) message("box_basic skipped (", m, "): ", conditionMessage(e))
  )
  tryCatch(
    ggsave(
      file.path(dir_peak, paste0("peak_", safe, ".png")),
      plot = plot_peak_timepoints(m),
      width = 10,
      height = 8,
      dpi = 300
    ),
    error = function(e) message("peak skipped (", m, "): ", conditionMessage(e))
  )
}

