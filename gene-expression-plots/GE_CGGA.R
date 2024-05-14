library(tidyverse)
library(matrixStats)
library(factoextra)
library(umap)
library(ggrepel)
library(pheatmap)

# ---- INPUTS ----
dir <- "~/path/to/project/"
setwd(dir)

groups_table <- "groups_v2.tsv"
counts_table <- "counts.tsv"

# ---- IMPORT TABLES ----

groups <- read_delim(groups_table) %>%
    # run this lines if you want to select samples coming from CGGA
    select(-GroupName) %>%
    filter(str_detect(Group, "CGGA"))

counts <- read_delim(counts_table) %>%
    select(all_of(c(groups$SampleName, "GENEID")))

print(str_glue("Number of unique GENE IDs: {length(unique(counts$GENEID))}
               Number of samples: {length(colnames(counts)) - 1}"))

# ---- FILTER BY VARIANCE ----
counts.matrix <- counts %>% 
    column_to_rownames("GENEID") %>%
    drop_na()%>%
    as.matrix()

variance <- rowVars(counts.matrix) / if_else(rowMeans(counts.matrix) != 0, rowMeans(counts.matrix), 1)

variance %>%
    as.data.frame() %>% 
    ggplot(aes(x = .)) +
    geom_histogram(binwidth = 3) +
    xlim(c(0, quantile(variance)[4]))

# take top 25% of genes in terms of variance
th <- quantile(variance)[3]
keep <- which(variance >= 15)

counts.matrix <- counts.matrix[keep, ]

# ---- PCA ----
PCs <- prcomp(t(counts.matrix), center = T, scale. = T)
idx <- match(colnames(counts.matrix), groups$SampleName)

fviz_pca_ind(
    PCs,
    col.ind = groups$Group[idx],
    geom = c("point"),
    pointsize = 4
) +
    scale_shape_manual(values = rep(20, nrow(PCs$x)))

# ---- UMAP ----
myumap <- umap(t(counts.matrix))

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

# ---- HEATMAP ----
group_annot <- groups %>%
    column_to_rownames("SampleName") %>%
    arrange(Group) %>%
    select(Group)

pheatmap(
    log2(counts.matrix + 1),  # para rotar, t() y cambiar annots (row por col y viceversa)
    show_rownames = F,
    show_colnames = F,
    cluster_cols = T, #si quieres clusterizar por columnas, es decir, muestras en este caso
    cluster_rows = F, #si quieres clusterizar por rows, es decir, evento en este caso
    clustering_distance_cols = "correlation",
    cellwidth = 2,
    annotation_col = group_annot
    # annotation_colors = list(Group = c(GSC = "red", hNSC = "blue"))
)
