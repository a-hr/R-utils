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
library(factoextra)


# ---- INPUTS ----
dir <- "~/path/to/alternative-splicing-plots"
groups_table <- "design.tab"
psi_table <- "PSI_TABLE.tab"

grupos_muestra <- c("BIOD_GSC") # muestras a mostrar en heatmaps/boxplots
# si de deja vaío, se usan todas

# ---- FORMATTING ----
setwd(dir)

groups <- read_delim(groups_table) %>%
    arrange(Group)

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

if (is.null(grupos_muestra)) {
    print("No specific groups selected, all will be used.")
    grupos_muestra <- unique(groups$GroupName)
}

# ---- K-means ----
# construct the matrix
samples <- groups %>%
    filter(Group %in% grupos_muestra) %>%
    filter(!str_detect(SampleName, "e")) %>%
    filter(!str_detect(SampleName, "GSC22_12_p6")) %>%
    pull(SampleName)

psi.matrix <- psi_table[, samples] %>% 
    scale() %>%
    as.data.frame() %>%
    drop_na() %>%
    t() %>%  # assuming columns => samples, rows => events
    as.matrix()

fviz_nbclust(psi.matrix, kmeans, method = "wss", k.max = 5)
fviz_nbclust(psi.matrix, kmeans, method = "silhouette", k.max = 5)

km <- kmeans(psi.matrix, centers = 2, nstart = 10)
fviz_cluster(km, data = psi.matrix, geom = "point") +
    theme_minimal()

clusters <- km$cluster

myumap <- umap(psi.matrix)

data.frame(
    x = myumap$layout[, 1],
    y = myumap$layout[, 2],
    k = clusters,
    samples = rownames(myumap$layout)
) %>%
    ggplot(aes(x=x, y=y, col = factor(clusters), label = samples)) +
    geom_point(size = 5) +
    geom_text_repel() +
    theme_void()
