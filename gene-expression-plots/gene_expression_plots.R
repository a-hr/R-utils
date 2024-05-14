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
dir <- "~/path/to/project/"
groups_table <- "GE_plots/groups_v2.tsv"
tpm_counts <- "GE_plots/TPM_COUNTS_AGGREGATED-lengthScaledTPM.tab"

grupos_muestra <- c("BIOD_BRAIN ", "BIOD_Astro", "CGGA_ASTR")  # muestras a mostrar en heatmaps/boxplots
gene_ratio <- 1  # mostrar el top N% de genes en funcion de la varianza (si es 1, se ignora el filtro con varianza)

# sacar los grupos de genes de una lista
#   Nota: si se deja la lista vacia, se usaran todos los genes (se pueden filtrar con gene_ratio)
target_genes <- c("MSI1", "MSI2", "SCYL3")

# sacar los grupos de genes de una tabla
target_genes <- read_delim("nombre-de-mi-tabla") %>%  # cambiar por el nombre del archivo con la lista de genes
    dplyr::select(GeneName) %>% # genes a mostrar (GeneName es el nombre de la columna donde estan)
    pull(GeneName)

# ---- FORMATTING ----
setwd(dir)
groups <- read_delim(groups_table) %>%
    dplyr::arrange(SampleName)

counts <- read_delim(tpm_counts) %>%
    dplyr::select(-GENEID)

if (is.null(target_genes)) {
    print("No target genes selected, all will be used")
    target_genes <- unique(counts$GENENAME)
}

if (is.null(grupos_muestra)) {
    print("No specific groups selected, all will be used.")
    grupos_muestra <- unique(groups$Group)
}



MY_THEME <-
    theme(
        text = element_text(family = "Roboto"),
        axis.text.x = element_text(angle = 35, vjust = .6),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_line(colour = "grey50"),
        panel.grid = element_line(color = "#b4aea9"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(linetype = "dashed"),
        panel.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
        plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
        legend.background = element_rect(fill = "#fbf9f4"),
        plot.title = element_text(
            family = "Roboto",
            size = 16,
            face = "bold",
            color = "#2a475e",
            margin = margin(b = 20)
        )
    )

# ---- BOXPLOT ----
# test
plot.data <- counts %>%
    dplyr::filter(GENENAME %in% target_genes) %>%
    pivot_longer(cols = which(colnames(.) != "GENENAME"), names_to = "SampleName", values_to = "Counts") %>%
    left_join(y = groups, by = "SampleName") %>%
    filter(Group %in% grupos_muestra)
comparisons <- combn(plot.data$Group %>% unique(), 2, simplify = FALSE)

# plot
plot.data %>% 
    ggplot(aes(x = Group, y = Counts + 1, col = Group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(size = .8, alpha = .8, position = position_jitterdodge()) +
    facet_wrap(~GENENAME, nrow = 2, ncol = 1, scales = "free_y") +
    guides(col = "none") +
    ylab("Counts") +
    scale_y_log10() +
    # comentar la linea stat_compare_means para quitar los t.tests
    # label "p.format" muestra p.val, "p.signif" muestra ns,*,**,***
    stat_compare_means(comparisons = comparisons, method = "t.test", label = "p.format") +
    MY_THEME

# ---- HEATMAP ----
muestras <- groups %>%
    filter(Group %in% grupos_muestra) %>%
    pull(SampleName) %>%
    unique()

hm.data <- counts %>%
    dplyr::filter(GENENAME %in% target_genes) %>%
    column_to_rownames("GENENAME") %>%
    select(all_of(muestras)) %>%
    as.matrix() %>% t()

col_annot <- groups %>%
    filter(SampleName %in% muestras) %>%
    column_to_rownames("SampleName") %>%
    arrange(Group)

g.order <- match(rownames(col_annot), rownames(hm.data))

hm.data <- hm.data[g.order,]
hm.log.data <- log2(hm.data + 1) %>% as.matrix()


# HEATMAP sin scaling, sin clustering
pheatmap(
    hm.data,
    show_rownames = F,
    show_colnames = T,
    cluster_cols = F,
    cluster_rows = F,
    cellwidth = 50,
    border_color = FALSE,
    annotation_row  = col_annot
)

# HEATMAP scaling por MUESTRA, sin clustering
pheatmap(
    hm.data,
    show_rownames = F,
    show_colnames = T,
    scale = "row", #asi me escala por row, es decir, en este caso por muestra. si pongo "column" en este caso escala por gen
    cluster_cols = F,
    cluster_rows = F,
    cellwidth = 50,
    border_color = FALSE,
    annotation_row  = col_annot
)

# HEATMAP scaling por GEN, sin clustering
pheatmap(
    hm.data,
    show_rownames = F,
    show_colnames = T,
    scale = "column", #asi me escala por gen, es decir, en este caso por muestra. si pongo "column" en este caso escala por gen
    cluster_cols = F,
    cluster_rows = F,
    cellwidth = 50,
    border_color = FALSE,
    annotation_row  = col_annot
)

# HEATMAP scaling LOG2, sin clustering
pheatmap(
    hm.log.data,
    show_rownames = F,
    show_colnames = T,
    scale = "none", 
    cluster_cols = F,
    cluster_rows = F,
    cellwidth = 50,
    border_color = FALSE,
    annotation_row  = col_annot
)

# HEATMAP scaling LOG2, con nombre por MUESTRA, sin clustering
pheatmap(
    hm.log.data,
    show_rownames = T,
    show_colnames = T,
    scale = "none", 
    cluster_cols = F,
    cluster_rows = F,
    cellwidth = 50,
    border_color = FALSE,
    annotation_row  = col_annot
)

# HEATMAP scaling LOG2, con clustering
pheatmap(
    hm.log.data,
    show_rownames = F,
    show_colnames = T,
    scale = "none", 
    cluster_cols = F, #si quieres clusterizar por columnas, es decir, gen en este caso
    cluster_rows = T, #si quieres clusterizar por rows, es decir, muestras en este caso
    clustering_distance_rows = "correlation",
    cellwidth = 50,
    border_color = FALSE,
    annotation_row  = col_annot
)

# ---- PCA ----

# construct the matrix
selected_cols <- groups %>%
    filter(Group %in% grupos_muestra) %>%
    pull(SampleName)

counts.matrix <- counts %>% 
    dplyr::filter(GENENAME %in% target_genes) %>%
    select(all_of(selected_cols)) %>% 
    t() %>%  # assuming columns => samples, rows => genes
    as.matrix()

# filter and scale the data
if (gene_ratio != 1) {
    variances <- colVars(counts.matrix) / colMaxs(counts.matrix)
    variances[is.na(variances)] <- 0
    threshold <- unname(quantile(abs(variances), probs = 1 - gene_ratio))
    keep <- which(variances >= threshold)  # keep genes with top variance
    
    counts.matrix <- log2(counts.matrix[, keep] + 1)
} else {
    counts.matrix <- log2(counts.matrix + 1)
}

# perform PCA
PCs <- prcomp(counts.matrix)

# variance explained by each PC (PC1 usually captures most of the noise-related data)
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
sample_groups <- sample_groups$Group

pca.plot1 <- data.frame(
    PC1 = PCs$x[, 1],
    PC2 = PCs$x[, 2],
    PC3 = PCs$x[, 3],
    Group = sample_groups
) %>%
    ggplot(aes(x = PC2, y = PC1, col = Group)) +
    geom_point() +
    theme_minimal()

pca.plot1
    
# ---- UMAP ----

# construct the matrix
selected_cols <- groups %>%
    filter(Group %in% grupos_muestra) %>%
    pull(SampleName)

counts.matrix <- counts %>% 
    dplyr::filter(GENENAME %in% target_genes) %>%
    select(all_of(selected_cols)) %>% 
    t() %>%  # assuming columns => samples, rows => genes
    as.matrix()

# scale the data
counts.matrix <- log2(counts.matrix + 1)

# perform dim reduction
myumap <- umap(counts.matrix, verbose=T)

# plot UMAP
sample_indices <- match(rownames(myumap$layout), groups$SampleName)

cluster_labels <- myumap$layout %>%
    as.data.frame() %>%
    rownames_to_column(var = "SampleName") %>%
    mutate(cluster = groups$Group[sample_indices]) %>%
    group_by(cluster) %>%
    summarize(x = mean(V1), y = mean(V2))

umap.df <- myumap$layout %>% 
    as.data.frame() %>%
    dplyr::rename("x" = "V1", "y" = "V2") %>%
    mutate(
        Group = groups$Group[sample_indices]
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

