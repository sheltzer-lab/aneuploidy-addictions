#!/usr/bin/env Rscript

library(ComplexHeatmap, quietly = TRUE)
library(reshape2, quietly = TRUE)
library(tidyverse, quietly = TRUE)
library(ggrepel, quietly = TRUE)

source(here::here("bin/functions.R"))
source(here::here("bin/constants.R"))

args <- commandArgs(trailingOnly = TRUE)

overallMutex <- lapply(unlist(stringr::str_split(args[6], ",")) %>% keep(~ str_detect(., "overall")), function(fh) {
    read.csv(fh, stringsAsFactors = FALSE)
}) %>%
    bind_rows() %>%
    mutate(arm = factor(arm,
        ordered = TRUE,
        levels = chroms
    ))

cancerMutexs <- lapply(unlist(stringr::str_split(args[6], ",")) %>% keep(~ !str_detect(., "overall")), function(fh) {
    read.csv(fh, stringsAsFactors = FALSE)
}) %>%
    bind_rows() %>%
    mutate(arm = factor(arm,
        ordered = TRUE,
        levels = chroms
    ))

mutations <- read.csv(args[3], stringsAsFactors = FALSE)
scores <- read.csv(args[4], stringsAsFactors = FALSE, check.names = FALSE) # check.names prevents naming as 'X1q' rather than '1q'
patients <- read.csv(args[5], stringsAsFactors = FALSE)

topN <- mutations %>%
    group_by(Hugo_Symbol) %>%
    summarise(n = n()) %>%
    slice_max(order_by = n, n = 25) %>%
    pluck("Hugo_Symbol")

z.mat <- overallMutex %>%
    filter(gene %in% topN) %>%
    acast(gene ~ arm, value.var = "Z")

png("overall_mutex_heatmap.png", width = 2560, height = 1361)
ht <- Heatmap(
    z.mat,
    name = "Z",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    heatmap_legend_param = list(
        legend_height = unit(1.5, "npc"),
        grid_width = unit(15, "mm"),
        title_gp = gpar(fontsize = 35, fontface = "bold"),
        labels_gp = gpar(fontsize = 35, fontface = "bold")
    ),
    rect_gp = gpar(col = "white", lwd = 2),
    row_names_gp = gpar(fontsize = 35, fontface = "bold"),
    column_names_gp = gpar(fontsize = 35, fontface = "bold")
)
draw(ht)
invisible(dev.off())

# For Fig1E show select overall examples of mutual exclusivity and co-occurrence
# OncoPrint plot
binary.mutated <- mutations %>%
    select(Hugo_Symbol, Tumor_Sample_Barcode) %>%
    acast(Hugo_Symbol ~ Tumor_Sample_Barcode, length) %>%
    {
        ifelse(. > 0, "MUT", NA)
    }
binary.gained <- scores %>%
    select(-PATIENT_ID) %>%
    mutate(across(where(is.numeric) &
        !c(SAMPLE_ID), function(x) {
        ifelse(x == 1, "GAIN", NA)
    })) %>%
    column_to_rownames("SAMPLE_ID") %>%
    as.matrix() %>%
    t()
matrix.intersection <- intersect(colnames(binary.mutated), colnames(binary.gained))
binary.together <- rbind(binary.mutated[, matrix.intersection], binary.gained[, matrix.intersection])
to.display <- c(
    "13q",
    "PTEN",
    "EGFR",
    "FOXA1",
    "MET",
    "TP53",
    "KRAS",
    "SMAD4",
    "NRAS"
)
fig1e.col <- c(COOC = "blue", MUTEX = "red", GAIN = "#FFC300")
fig1e.alter_fun <- list(
    background = alter_graphic("rect", fill = "#ECECEC"),
    COOC = alter_graphic("rect", fill = fig1e.col["COOC"]),
    MUTEX = alter_graphic("rect", fill = fig1e.col["MUTEX"]),
    GAIN = alter_graphic("rect", fill = fig1e.col["GAIN"])
)
m <- binary.together[to.display, intersect(colnames(binary.together), patients[, "SAMPLE_ID"])] %>%
    t() %>%
    as.data.frame() %>%
    mutate(across(c("PTEN", "EGFR", "FOXA1", "MET"), ~ str_replace_all(., "MUT", "MUTEX"))) %>%
    mutate(across(c("TP53", "KRAS", "NRAS", "SMAD4"), ~ str_replace_all(., "MUT", "COOC"))) %>%
    as.matrix() %>%
    t()
png("overall_examples.png", width = 2560, height = 1361)
op <- oncoPrint(
    m,
    alter_fun = fig1e.alter_fun,
    col = fig1e.col,
    show_heatmap_legend = FALSE,
    remove_empty_columns = TRUE,
    # remove_empty_rows = TRUE,
    row_order = to.display,
    row_split = factor(c("Arm", rep("Mutually Exclusive", 4), rep("Co-occurrence", 4)), levels = c("Arm", "Mutually Exclusive", "Co-occurrence")),
    cluster_row_slices = FALSE,
    cluster_column_slices = FALSE,
    top_annotation = NULL,
    right_annotation = NULL,
    row_gap = unit(5, "mm"),
    row_names_gp = gpar(fontsize = 50, family = "bold"),
    row_title_gp = gpar(fontsize = 55, family = "bold"),
    pct_gp = gpar(fontsize = 40),
    show_pct = FALSE,
    row_title = NULL
)
draw(op)
invisible(dev.off())


# Begin building Supplemental Figure 1
# Cancer-specific heatmaps and oncoPrints to match Fig 1D and Fig 1E
makeHeatmap <- function(z.mat, ...) {
    Heatmap(
        z.mat,
        name = "Z",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        heatmap_legend_param = list(
            legend_height = unit(0.25, "npc"),
            grid_width = unit(3, "mm"),
            title_gp = gpar(fontsize = 6, fontface = "bold"),
            labels_gp = gpar(fontsize = 6, fontface = "bold"),
            title = ""
        ),
        rect_gp = gpar(col = "white", lwd = 1),
        row_names_gp = gpar(fontsize = 6, fontface = "bold"),
        column_names_gp = gpar(fontsize = 6, fontface = "bold"),
        col = circlize::colorRamp2(c(-20, 0, 20), c("blue", "white", "red")),
        ...
    )
}
makeOncoPrint <- function(m, ...) {
    oncoPrint(
        m,
        alter_fun = list(
            background = alter_graphic("rect", fill = "#ECECEC"),
            COOC = alter_graphic("rect", fill = "blue"),
            MUTEX = alter_graphic("rect", fill = "red"),
            GAIN = alter_graphic("rect", fill = "#FFC300")
        ),
        show_heatmap_legend = FALSE,
        remove_empty_columns = TRUE,
        # remove_empty_rows = TRUE,
        row_order = to.display,
        row_split = factor(
            c("Arm", rep("Mutually Exclusive", 4), rep("Co-occurrence", 4)),
            levels = c("Arm", "Mutually Exclusive", "Co-occurrence")
        ),
        cluster_row_slices = FALSE,
        cluster_column_slices = FALSE,
        top_annotation = NULL,
        right_annotation = NULL,
        row_gap = unit(1, "mm"),
        row_names_gp = gpar(fontsize = 6, family = "bold"),
        row_title_gp = gpar(fontsize = 6, family = "bold"),
        pct_gp = gpar(fontsize = 3),
        show_pct = FALSE,
        row_title = NULL,
        ...
    )
}

pancreatic.heatmap <- makeHeatmap(
    cancerMutexs %>%
        filter(gene %in% topN, cancer == "pancreatic_cancer") %>%
        acast(gene ~ arm, value.var = "Z")
)
lungcancer.heatmap <- makeHeatmap(
    cancerMutexs %>%
        filter(gene %in% topN, cancer == "non_small_cell_lung_cancer") %>%
        acast(gene ~ arm, value.var = "Z")
)
colorectal.heatmap <- makeHeatmap(
    cancerMutexs %>%
        filter(gene %in% topN, cancer == "colorectal_cancer") %>%
        acast(gene ~ arm, value.var = "Z")
)
breast.heatmap <- makeHeatmap(
    cancerMutexs %>%
        filter(gene %in% topN, cancer == "breast_cancer") %>%
        acast(gene ~ arm, value.var = "Z")
)

selections <- tribble(
    ~cancer.title, ~cancer.clean, ~select.arm,
    "Pancreatic Cancer", "pancreatic_cancer", "18q",
    "Non-Small Cell Lung Cancer", "non_small_cell_lung_cancer", "7p",
    "Colorectal Cancer", "colorectal_cancer", "13q",
    "Breast Cancer", "breast_cancer", "5q"
) %>%
    rowwise() %>%
    mutate(
        genes.cooccurring = cancerMutexs %>% filter(cancer == cancer.clean, arm == select.arm) %>% arrange(Z) %>% pluck("gene") %>% head(n = 4) %>% list(),
        genes.mutex = cancerMutexs %>% filter(cancer == cancer.clean, arm == select.arm) %>% arrange(desc(Z)) %>% pluck("gene") %>% head(n = 4) %>% list()
    )

# "pancreatic_cancer" oncoprint
subset.patients <- patients[which(patients$CANCER_TYPE == "Pancreatic Cancer"), "SAMPLE_ID"]
binary.mutated <- mutations %>%
    filter(Tumor_Sample_Barcode %in% subset.patients) %>%
    select(Hugo_Symbol, Tumor_Sample_Barcode) %>%
    acast(Hugo_Symbol ~ Tumor_Sample_Barcode, length) %>%
    {
        ifelse(. > 0, "MUT", NA)
    }
binary.gained <- scores %>%
    select(-PATIENT_ID) %>%
    filter(SAMPLE_ID %in% subset.patients) %>%
    mutate(across(where(is.numeric) &
        !c(SAMPLE_ID), function(x) {
        ifelse(x == 1, "GAIN", NA)
    })) %>%
    column_to_rownames("SAMPLE_ID") %>%
    as.matrix() %>%
    t()
matrix.intersection <- intersect(colnames(binary.mutated), colnames(binary.gained))
binary.together <- rbind(binary.mutated[, matrix.intersection], binary.gained[, matrix.intersection])
arm.selection <- selections %>%
    filter(cancer.clean == "pancreatic_cancer") %>%
    pluck("select.arm")
mutex.selection <- selections %>%
    filter(cancer.clean == "pancreatic_cancer") %>%
    pluck("genes.mutex") %>%
    unlist()
coocurrence.selection <- selections %>%
    filter(cancer.clean == "pancreatic_cancer") %>%
    pluck("genes.cooccurring") %>%
    unlist()
to.display <- c(arm.selection, mutex.selection, coocurrence.selection)
m <- binary.together[to.display, intersect(colnames(binary.together), patients[, "SAMPLE_ID"])] %>%
    t() %>%
    as.data.frame() %>%
    mutate(across(all_of(mutex.selection), ~ str_replace_all(., "MUT", "MUTEX"))) %>%
    mutate(across(all_of(coocurrence.selection), ~ str_replace_all(., "MUT", "COOC"))) %>%
    as.matrix() %>%
    t()
pancreatic.onoprint <- makeOncoPrint(m)

# "non_small_cell_lung_cancer" oncoprint
subset.patients <- patients[which(patients$CANCER_TYPE == "Non-Small Cell Lung Cancer"), "SAMPLE_ID"]
binary.mutated <- mutations %>%
    filter(Tumor_Sample_Barcode %in% subset.patients) %>%
    select(Hugo_Symbol, Tumor_Sample_Barcode) %>%
    acast(Hugo_Symbol ~ Tumor_Sample_Barcode, length) %>%
    {
        ifelse(. > 0, "MUT", NA)
    }
binary.gained <- scores %>%
    select(-PATIENT_ID) %>%
    filter(SAMPLE_ID %in% subset.patients) %>%
    mutate(across(where(is.numeric) &
        !c(SAMPLE_ID), function(x) {
        ifelse(x == 1, "GAIN", NA)
    })) %>%
    column_to_rownames("SAMPLE_ID") %>%
    as.matrix() %>%
    t()
matrix.intersection <- intersect(colnames(binary.mutated), colnames(binary.gained))
binary.together <- rbind(binary.mutated[, matrix.intersection], binary.gained[, matrix.intersection])
arm.selection <- selections %>%
    filter(cancer.clean == "non_small_cell_lung_cancer") %>%
    pluck("select.arm")
mutex.selection <- selections %>%
    filter(cancer.clean == "non_small_cell_lung_cancer") %>%
    pluck("genes.mutex") %>%
    unlist()
coocurrence.selection <- selections %>%
    filter(cancer.clean == "non_small_cell_lung_cancer") %>%
    pluck("genes.cooccurring") %>%
    unlist()
to.display <- c(arm.selection, mutex.selection, coocurrence.selection)
m <- binary.together[to.display, intersect(colnames(binary.together), patients[, "SAMPLE_ID"])] %>%
    t() %>%
    as.data.frame() %>%
    mutate(across(all_of(mutex.selection), ~ str_replace_all(., "MUT", "MUTEX"))) %>%
    mutate(across(all_of(coocurrence.selection), ~ str_replace_all(., "MUT", "COOC"))) %>%
    as.matrix() %>%
    t()
lungcancer.onoprint <- makeOncoPrint(m)


# "colorectal_cancer" oncoprint
subset.patients <- patients[which(patients$CANCER_TYPE == "Colorectal Cancer"), "SAMPLE_ID"]
binary.mutated <- mutations %>%
    filter(Tumor_Sample_Barcode %in% subset.patients) %>%
    select(Hugo_Symbol, Tumor_Sample_Barcode) %>%
    acast(Hugo_Symbol ~ Tumor_Sample_Barcode, length) %>%
    {
        ifelse(. > 0, "MUT", NA)
    }
binary.gained <- scores %>%
    select(-PATIENT_ID) %>%
    filter(SAMPLE_ID %in% subset.patients) %>%
    mutate(across(where(is.numeric) &
        !c(SAMPLE_ID), function(x) {
        ifelse(x == 1, "GAIN", NA)
    })) %>%
    column_to_rownames("SAMPLE_ID") %>%
    as.matrix() %>%
    t()
matrix.intersection <- intersect(colnames(binary.mutated), colnames(binary.gained))
binary.together <- rbind(binary.mutated[, matrix.intersection], binary.gained[, matrix.intersection])
arm.selection <- selections %>%
    filter(cancer.clean == "colorectal_cancer") %>%
    pluck("select.arm")
mutex.selection <- selections %>%
    filter(cancer.clean == "colorectal_cancer") %>%
    pluck("genes.mutex") %>%
    unlist()
coocurrence.selection <- selections %>%
    filter(cancer.clean == "colorectal_cancer") %>%
    pluck("genes.cooccurring") %>%
    unlist()
to.display <- c(arm.selection, mutex.selection, coocurrence.selection)
m <- binary.together[to.display, intersect(colnames(binary.together), patients[, "SAMPLE_ID"])] %>%
    t() %>%
    as.data.frame() %>%
    mutate(across(all_of(mutex.selection), ~ str_replace_all(., "MUT", "MUTEX"))) %>%
    mutate(across(all_of(coocurrence.selection), ~ str_replace_all(., "MUT", "COOC"))) %>%
    as.matrix() %>%
    t()
colorectal.onoprint <- makeOncoPrint(m)

# "breast_cancer" oncoprint
subset.patients <- patients[which(patients$CANCER_TYPE == "Breast Cancer"), "SAMPLE_ID"]
binary.mutated <- mutations %>%
    filter(Tumor_Sample_Barcode %in% subset.patients) %>%
    select(Hugo_Symbol, Tumor_Sample_Barcode) %>%
    acast(Hugo_Symbol ~ Tumor_Sample_Barcode, length) %>%
    {
        ifelse(. > 0, "MUT", NA)
    }
binary.gained <- scores %>%
    select(-PATIENT_ID) %>%
    filter(SAMPLE_ID %in% subset.patients) %>%
    mutate(across(where(is.numeric) &
        !c(SAMPLE_ID), function(x) {
        ifelse(x == 1, "GAIN", NA)
    })) %>%
    column_to_rownames("SAMPLE_ID") %>%
    as.matrix() %>%
    t()
matrix.intersection <- intersect(colnames(binary.mutated), colnames(binary.gained))
binary.together <- rbind(binary.mutated[, matrix.intersection], binary.gained[, matrix.intersection])
arm.selection <- selections %>%
    filter(cancer.clean == "breast_cancer") %>%
    pluck("select.arm")
mutex.selection <- selections %>%
    filter(cancer.clean == "breast_cancer") %>%
    pluck("genes.mutex") %>%
    unlist()
coocurrence.selection <- selections %>%
    filter(cancer.clean == "breast_cancer") %>%
    pluck("genes.cooccurring") %>%
    unlist()
to.display <- c(arm.selection, mutex.selection, coocurrence.selection)
m <- binary.together[to.display, intersect(colnames(binary.together), patients[, "SAMPLE_ID"])] %>%
    t() %>%
    as.data.frame() %>%
    mutate(across(all_of(mutex.selection), ~ str_replace_all(., "MUT", "MUTEX"))) %>%
    mutate(across(all_of(coocurrence.selection), ~ str_replace_all(., "MUT", "COOC"))) %>%
    as.matrix() %>%
    t()
breast.onoprint <- makeOncoPrint(m)

plt <- cowplot::plot_grid(
    cowplot::plot_grid(
        grid.grabExpr(draw(pancreatic.heatmap)),
        grid.grabExpr(draw(lungcancer.heatmap)),
        grid.grabExpr(draw(colorectal.heatmap)),
        grid.grabExpr(draw(breast.heatmap)),
        ncol = 1
    ),
    NULL,
    cowplot::plot_grid(
        grid.grabExpr(draw(pancreatic.onoprint)),
        NULL,
        grid.grabExpr(draw(lungcancer.onoprint)),
        NULL,
        grid.grabExpr(draw(colorectal.onoprint)),
        NULL,
        grid.grabExpr(draw(breast.onoprint)),
        NULL,
        ncol = 1,
        rel_heights = c(rep(c(1, 0.065), 4))
    ),
    ncol = 3,
    rel_widths = c(1.0, 0.15, 1.0)
)
cowplot::save_plot(filename = "supplemental-figure.png", plot = plt, base_height = 11.7, base_width = 8.3)
