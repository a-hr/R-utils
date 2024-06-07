library(maser)
library(tidyverse)
library(factoextra)
library(patchwork)


# ---- Read Data ----
dir <- "~/Library/Mobile Documents/com~apple~CloudDocs/Biodonostia/flecanda"
setwd(dir)

dir.create("tables")
dir.create("plots")
dir.create("plots/volcanos")
dir.create("plots/PCAs")

# load design dataframe
experiments <- list.files("rmats/")

# load rMATS data
masers <- list()

for (exp in experiments) {
    conditions <- exp %>% str_remove("rmats_") %>% str_split("_")
    conditions <- conditions[[1]]
    events <- maser(path = str_glue("rmats/{exp}/rmats_output"), cond_labels = conditions, ftype = c("JC"))
    masers[[exp]] <- filterByCoverage(events, avg_reads = 20)
}

MY_THEME <-
    theme(
        text = element_text(family = "Roboto"),
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

# ---- Volcano ----
events <- c("A3SS", "A5SS", "SE", "RI", "MXE")

for (exp in experiments) {
    conditions <- exp %>% str_remove("rmats_") %>% str_split("_")
    conditions <- conditions[[1]]
    volcanos <-
        lapply(events, function(t) {
            volcano(masers[[exp]], type = t, deltaPSI = .1, fdr = .05) + ggtitle(t)
        })
    
    y.max <- max(sapply(volcanos, function(p) max(p$data$log10pval)))
    y.min <- max(sapply(volcanos, function(p) min(p$data$log10pval)))
    
    legend.levels <- levels(lapply(volcanos, function(p) p$data$Status)[[1]])
    color_scale <- c("gray", "tomato2", "navy")
    names(color_scale) <- legend.levels
    
    v.n <- wrap_plots(volcanos, guides = "collect") +
        guide_area() +
        plot_annotation(title = str_glue("{conditions[1]} vs {conditions[2]}")) &
        scale_color_manual(values = color_scale) &
        ylim(c(y.min, y.max)) &
        MY_THEME
    
    v.n
    
    ggsave(
        filename = str_glue("plots/volcanos/volcano_{conditions[1]}_vs_{conditions[2]}.png"),
        plot = v.n,
        device = "png",
        width = 400,
        height = 250,
        units = "mm",
        dpi = 320,
        bg = "white"
    )
    
}

# ---- PCA ----
events <- c("A3SS", "A5SS", "SE", "RI", "MXE")

for (exp in experiments) {
    conditions <- exp %>% str_remove("rmats_") %>% str_split("_")
    conditions <- conditions[[1]]
    
    psis <- do.call(
        rbind,
        lapply(events, function(t)
            PSI(masers[[exp]], type = t))
        ) %>%
        as.data.frame() %>%
        drop_na() %>%
        scale() %>%
        t() %>%
        as.matrix()
    
    # compute PCA
    PCs <- prcomp(psis, center = T, scale. = F)
    
    # visualize components
    # pca.1 <- fviz_eig(PCs)
    
    # visualize PCA
    pca <- fviz_pca_ind(
        PCs,
        col.ind = c(rep(masers[[exp]]@conditions[1], 3), rep(masers[[exp]]@conditions[2], 3)),
        geom = c("point"),
        pointsize = 4
    ) +
        scale_shape_manual(values = rep(20, nrow(psis))) +
        ggtitle(str_glue("{conditions[1]} vs {conditions[2]}")) +
        MY_THEME
    
    ggsave(
        filename = str_glue("plots/PCAs/PCA_{conditions[1]}_vs_{conditions[2]}.png"),
        plot = pca,
        device = "png",
        width = 250,
        height = 250,
        units = "mm",
        dpi = 320,
        bg = "white"
    )
}

# ---- Tables ----
tables <- list()

# All events per sample (only common columns included)
for (exp.maser in masers) {
    all.events <- lapply(as.list(events), function(event)
        summary(exp.maser, type = event) %>% select(all_of(
            c(
                "ID",
                "GeneID",
                "geneSymbol",
                "FDR",
                "IncLevelDifference",
                "PSI_1",
                "PSI_2",
                "Chr",
                "Strand"
            )
        )) %>%
            filter(FDR <= .05, abs(IncLevelDifference) >= .1) %>%
            mutate(Type = event)
    )
    contrast <- paste(exp.maser@conditions, collapse = "_vs_")
    tables[[contrast]] <- do.call(rbind, all.events)
}

for (ii in seq_along(tables)) {
    table.name <- names(tables)[ii]
    write_delim(
        tables[[ii]],
        file = str_glue("tables/ALL_SIGNIF_EVENTS_{table.name}.tab"),
        delim = "\t"
    )
}

# One table per event and sample (full tables)
tables <- list()
for (exp.maser in masers) {
    event.tables <- lapply(as.list(events), function(event)
        summary(exp.maser, type = event) %>%
            filter(FDR <= .05, abs(IncLevelDifference) >= .1)
    )
    names(event.tables) <- events
    
    contrast <- paste(exp.maser@conditions, collapse = "_vs_")
    tables[[contrast]] <- event.tables
}

for (ii in seq_along(tables)) {
    exp.name <- names(tables)[ii]
    dir.create(str_glue("tables/{exp.name}"))
    for (jj in seq_along(tables[[ii]])) {
        event.name <- names(tables[[ii]])[jj]
        write_delim(
            tables[[ii]][[jj]],
            file = str_glue("tables/{exp.name}/{exp.name}_{event.name}.tab"),
            delim = "\t"
        )   
    }
}
