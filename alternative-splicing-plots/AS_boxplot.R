library(dplyr)
library(tidyverse)
library(tximport)
library(DESeq2)
library(umap)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(RColorBrewer)
library(pheatmap)
library(ggrepel)
library(patchwork)

# ---- INPUTS ----
dir <- "~/path/to/alternative-splicing-plots"
setwd(dir)

groups_table <- "groups_v2.csv"
psi_table <- "PSI_TABLE.tab"

grupos_muestra <- c("BIOD_GB", "BIOD_BRAIN ") # muestras a mostrar en heatmaps/boxplots
desired_order_of_groups <- c("BIOD_BRAIN ", "BIOD_GB")
group_colors <-c(BIOD_GB ="blue", `BIOD_BRAIN ` = "deeppink")

# sacar los eventos de una lista # Nota: si se deja la lista vacia, se usaran todos los eventos
selected_events <- c("HsaEX0047869", "HsaEX0020272")

# ---- FORMATTING ----
groups <- read_delim(groups_table) %>%
  dplyr::arrange(SampleName)

psi_table <- read_delim(
  psi_table,
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE,
  show_col_types = FALSE
) %>%
    dplyr::filter(EVENT %in% selected_events) %>%
mutate(
    TYPE = case_when(
      grepl("EX", EVENT) ~ "EX",
      grepl("INT", EVENT) ~ "INT",
      grepl("ALT", EVENT) ~ "ALT"
    )
  ) %>%
  drop_na()

# ---- BOXPLOT, grupos ----
## Get samples
samples <- groups %>%
  filter(Group %in% grupos_muestra) %>%
  pull(SampleName)

# test
plot.data <- psi_table[, c("EVENT", samples)] %>% 
    pivot_longer(cols = which(colnames(.) != "EVENT"), names_to = "SampleName", values_to = "PSI") %>%
    left_join(x = groups, by = "SampleName") %>%
    filter(Group %in% grupos_muestra)

comparisons <- combn(plot.data$Group %>% unique(), 2, simplify = FALSE)

# Ensure 'Group' is a factor with desired order
plot.data$Group <- factor(plot.data$Group, levels = desired_order_of_groups)

# plot
plot.data %>%
  ggplot(aes(x = Group, y = PSI + 1, col = Group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = .8, alpha = .8, position = position_jitterdodge()) +
  facet_wrap(~EVENT, nrow = 2, ncol = 4, scales = "free_y") +
  guides(col = "none") +
  ylab("Counts") +
  scale_y_log10() +
  # comentar la linea stat_compare_means para quitar los t.tests
  # label "p.format" muestra p.val, "p.signif" muestra ns,*,**,***
  stat_compare_means (comparisons = comparisons, method = "t.test", label = "p.signif") +
  scale_color_manual(values = group_colors)# Assign colors to groups
