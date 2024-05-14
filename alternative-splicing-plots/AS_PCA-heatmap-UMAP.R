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
dir <- "~/Library/Mobile Documents/com~apple~CloudDocs/Biodonostia/GB/AS_plots"
groups_table <- "design.tab"
psi_table <- "PSI_TABLE.tab"

grupos_muestra <- c("GSC", "hNSC") # muestras a mostrar en heatmaps/boxplots
# si de deja vacío, se usan todas

# sacar los eventos de una lista 
# Nota: si se deja la lista vacia, se usaran todos los eventos
selected_events <- c("MSI1", "MSI2", "SCYL3")

# sacar los eventos de una tabla
selected_events <- read_delim(paste0(dir, "/", "200411_GSCvshNSC_PCT50_DiffEXpsi_table_only_up.csv")) %>%  # cambiar por el nombre del archivo con la lista de eventos
    dplyr::select(EVENT) %>% # genes a mostrar (GeneName es el nombre de la columna donde estan)
    pull(EVENT)

# ---- FORMATTING ----
setwd(dir)

groups <- read_delim(groups_table) %>%
    arrange(GroupName)

psi_table <- read_delim(
    psi_table,
    delim = "\t",
    escape_double = FALSE,
    trim_ws = TRUE,
    show_col_types = FALSE
) %>%
    # dplyr::filter(EVENT %in% selected_events) %>%
    mutate(
        TYPE = case_when(
            grepl("EX", EVENT) ~ "EX",
            grepl("INT", EVENT) ~ "INT",
            grepl("ALT", EVENT) ~ "ALT"
        )
    ) %>%
    column_to_rownames('EVENT')

if (is.null(selected_events)) {
    print("No target events selected, all will be used")
    selected_events <- rownames(psi_table)
}

if (is.null(grupos_muestra)) {
    print("No specific groups selected, all will be used.")
    grupos_muestra <- unique(groups$GroupName)
}

# ---- PCA ----
# construct the matrix
samples <- groups %>%
    filter(GroupName %in% grupos_muestra) %>%
    pull(SampleName)

psi.matrix <- psi_table[selected_events, samples] %>% 
    t() %>%  # assuming columns => samples, rows => events
    as.matrix() %>%
    scale()

# perform PCA
PCs <- prcomp(psi.matrix)

# variance explained by each PC (PC1 usually captures most of the noise-related variance)
data.frame(
    variance = PCs$sdev/sum(PCs$sdev) * 100,
    PC = paste0("PC", 1:length(PCs$sdev))
) %>%
    dplyr::filter(row_number() <= 10) %>%
    mutate(PC = factor(PC, levels = PC)) %>%
    ggplot(aes(x = PC, y = variance)) +
    geom_bar(stat = "identity", fill = "skyblue") +
    theme_minimal()

# PCA plots
sample_groups <- merge(
    as.data.frame(rownames(PCs$x)), 
    groups, 
    by.x = "rownames(PCs$x)",
    by.y = "SampleName"
)
sample_groups <- sample_groups$GroupName

pca.plot <- data.frame(
    PC1 = PCs$x[, 1],
    PC2 = PCs$x[, 2],
    PC3 = PCs$x[, 3],
    Group = sample_groups
) %>%
    ggplot(aes(x = PC1, y = PC2, col = Group)) +
    geom_point() +
    theme_minimal()

pca.plot

# ---- HEATMAP ----
## Get samples
samples <- groups %>%
    filter(GroupName %in% grupos_muestra) %>%
    pull(SampleName)

## Get col and row annotations
event_annot <- psi_table %>%
    select(TYPE)

psi_table <- psi_table %>%
    dplyr::select(-any_of(c("EVENT", "TYPE")))

psi_table <- psi_table[, samples] # reorder table to match annotations

group_annot <- groups %>%
    filter(SampleName %in% samples) %>%
    column_to_rownames("SampleName") %>%
    arrange(GroupName) %>%
    dplyr::rename("Group" = "GroupName") %>%
    select(Group)

# HEATMAP sin scaling, sin clustering
pheatmap(
    psi_table,  # para rotar, t() y cambiar annots (row por col y viceversa)
    show_rownames = F,
    show_colnames = T,
    scale = "row",
    cluster_cols = F, #si quieres clusterizar por columnas, es decir, muestras en este caso
    cluster_rows = F, #si quieres clusterizar por rows, es decir, evento en este caso
    clustering_distance_cols = "correlation",
    cellwidth = 10,
    border_color = FALSE,
    annotation_row = event_annot,
    annotation_col = group_annot,
    annotation_colors = list(Group = c(GSC = "red", hNSC = "blue"))
)

# ---- UMAP ----
# construct the matrix
samples <- groups %>%
    filter(GroupName %in% grupos_muestra) %>%
    pull(SampleName)

psi.matrix <- psi_table[, samples] %>% 
    t() %>%  # assuming columns => samples, rows => events
    as.matrix()

# perform dim reduction
myumap <- umap(psi.matrix, verbose=T)

# plot UMAP
sample_indices <- match(rownames(myumap$layout), groups$SampleName)

cluster_labels <- myumap$layout %>%
    as.data.frame() %>%
    rownames_to_column(var = "SampleName") %>%
    mutate(cluster = groups$GroupName[sample_indices]) %>%
    group_by(cluster) %>%
    summarize(x = mean(V1), y = mean(V2))

umap.df <- myumap$layout %>% 
    as.data.frame() %>%
    dplyr::rename("x" = "V1", "y" = "V2") %>%
    mutate(
        Group = groups$GroupName[sample_indices]
    )

ggplot() +
    geom_point(data = umap.df, aes(x = x, y = y, col = Group)) +
    # comentar geom_label_repel para quitar los títulos
    geom_label_repel(
        data = cluster_labels, 
        aes(x = x, y = y, label = cluster, col = cluster), 
        alpha = .8
    ) +
    theme_void()
