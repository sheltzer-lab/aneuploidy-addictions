#!/usr/bin/env Rscript

library(assertr, quietly = TRUE)
library(tidyverse, quietly = TRUE)

# Determine the path to local shared R definitions
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)

source(file.path(script.basename, "constants.R"))
source(file.path(script.basename, "functions.R"))

args <- commandArgs(trailingOnly = TRUE)

study <- args[1]

patientsCutoffCna <- as.numeric(args[4])
patientsCutoffCnaRelative <- as.numeric(args[6])

samples <- read.csv(args[5], stringsAsFactors = FALSE) %>%
    verify(has_all_names("PATIENT_ID", "STATUS", "MONTHS", "SAMPLE_ID", "CANCER_TYPE")) %>%
    select(PATIENT_ID, SAMPLE_ID, CANCER_TYPE) %>%
    assert(is_uniq, PATIENT_ID, SAMPLE_ID)

mutation <- read.csv(args[2], stringsAsFactors = FALSE) %>%
    verify(has_all_names("Hugo_Symbol", "Tumor_Sample_Barcode", "Variant_Classification")) %>%
    count(Tumor_Sample_Barcode, Hugo_Symbol)

score <- read_csv(args[3], col_types = cols(.default = "c")) %>%
    verify(has_all_names("PATIENT_ID", "SAMPLE_ID")) %>%
    assert(is_uniq, PATIENT_ID, SAMPLE_ID) %>%
    select(-PATIENT_ID) %>%
    mutate(across(!c(SAMPLE_ID), as.numeric)) %>%
    assert(within_bounds(-1, 1), any_of(chroms))

# Produce mutual exclusivity tables by cancer
for (cancer in unique(samples$CANCER_TYPE)) {
    cleanCancerStr <- str_replace_all(cancer, "[^a-zA-Z0-9]+", "_") %>% str_to_lower()
    s <- samples %>%
        filter(CANCER_TYPE == cancer) %>%
        pluck("SAMPLE_ID") %>%
        unique()
    cutoff <- min(
        patientsCutoffCna,
        length(s) * patientsCutoffCnaRelative
    )

    # There are enough samples to surpass the cutoff
    if (cutoff <= length(s)) {
        subMutations <- mutation %>%
            filter(Tumor_Sample_Barcode %in% s) %>%
            pivot_wider(
                id_cols = Tumor_Sample_Barcode,
                names_from = Hugo_Symbol,
                values_from = n,
                values_fill = 0,
                values_fn = function(x) {
                    x >= 1
                }
            ) %>%
            column_to_rownames("Tumor_Sample_Barcode") %>%
            as.matrix()

        subScores <- score %>%
            filter(SAMPLE_ID %in% s)

        gain <- subScores %>%
            mutate(across(where(is.numeric) &
                !c(SAMPLE_ID), function(x) {
                x == 1
            })) %>%
            column_to_rownames("SAMPLE_ID") %>%
            as.matrix()

        gain.matrices <- mutexclusivity(subMutations, gain, minimumPatientsWcna = cutoff)
        if (dim(gain.matrices$neither)[[2]] != 0) {
            expand.grid(
                arm = colnames(gain.matrices$neither),
                gene = rownames(gain.matrices$neither)
            ) %>%
                as_tibble() %>%
                rowwise() %>%
                mutate(
                    neither = gain.matrices$neither[gene, arm],
                    both = gain.matrices$both[gene, arm],
                    gene.not.cna = gain.matrices$gene.not.cna[gene, arm],
                    cna.not.gene = gain.matrices$cna.not.gene[gene, arm],
                    pval = gain.matrices$fisher.p[gene, arm],
                    odds.ratio = gain.matrices$fisher.odds[gene, arm]
                ) %>%
                ungroup() %>%
                mutate(
                    study = study,
                    cancer = cleanCancerStr,
                    neg.log2.odds = -log2(odds.ratio)
                ) %>%
                mutate(
                    Z = sign(neg.log2.odds) * p.to.Z(pval) # Z of p-value with direction sign
                ) %>%
                group_by(study, cancer, arm) %>%
                mutate(
                    padj = p.adjust(pval, method = "fdr"),
                ) %>%
                ungroup() %>%
                select(study, cancer, arm, gene, neither, both, gene.not.cna, cna.not.gene, pval, padj, odds.ratio, neg.log2.odds, Z) %>%
                arrange(padj) %>%
                verify(has_all_names(
                    "study",
                    "cancer",
                    "arm",
                    "gene",
                    "neither",
                    "both",
                    "gene.not.cna",
                    "cna.not.gene",
                    "pval",
                    "padj",
                    "odds.ratio",
                    "neg.log2.odds",
                    "Z"
                )) %>%
                write.csv(paste0(study, "--", cleanCancerStr, "--arm-gain-mutexclusivity.csv"),
                    row.names = FALSE,
                    na = "",
                    quote = TRUE
                )
        }
    }
}


# Produce pancancer mutual exclusivity table
s <- samples %>%
    pluck("SAMPLE_ID") %>%
    unique()
cutoff <- min(
    patientsCutoffCna,
    length(s) * patientsCutoffCnaRelative
)
# There are enough samples to surpass the cutoff
if (cutoff <= length(s)) {
    mutMatrix <- mutation %>%
        pivot_wider(
            id_cols = Tumor_Sample_Barcode,
            names_from = Hugo_Symbol,
            values_from = n,
            values_fill = 0,
            values_fn = function(x) {
                x >= 1
            }
        ) %>%
        column_to_rownames("Tumor_Sample_Barcode") %>%
        as.matrix()

    gain <- score %>%
        mutate(across(where(is.numeric) &
            !c(SAMPLE_ID), function(x) {
            x == 1
        })) %>%
        column_to_rownames("SAMPLE_ID") %>%
        as.matrix()

    gain.matrices <- mutexclusivity(mutMatrix, gain, minimumPatientsWcna = cutoff)
    if (dim(gain.matrices$neither)[[2]] != 0) {
        expand.grid(
            arm = colnames(gain.matrices$neither),
            gene = rownames(gain.matrices$neither)
        ) %>%
            as_tibble() %>%
            rowwise() %>%
            mutate(
                neither = gain.matrices$neither[gene, arm],
                both = gain.matrices$both[gene, arm],
                gene.not.cna = gain.matrices$gene.not.cna[gene, arm],
                cna.not.gene = gain.matrices$cna.not.gene[gene, arm],
                pval = gain.matrices$fisher.p[gene, arm],
                odds.ratio = gain.matrices$fisher.odds[gene, arm]
            ) %>%
            ungroup() %>%
            mutate(
                study = study,
                cancer = study,
                neg.log2.odds = -log2(odds.ratio)
            ) %>%
            mutate(
                Z = sign(neg.log2.odds) * p.to.Z(pval) # Z of p-value with direction sign
            ) %>%
            group_by(study, cancer, arm) %>%
            mutate(
                padj = p.adjust(pval, method = "fdr"),
            ) %>%
            ungroup() %>%
            select(study, cancer, arm, gene, neither, both, gene.not.cna, cna.not.gene, pval, padj, odds.ratio, neg.log2.odds, Z) %>%
            arrange(padj) %>%
            verify(has_all_names(
                "study",
                "cancer",
                "arm",
                "gene",
                "neither",
                "both",
                "gene.not.cna",
                "cna.not.gene",
                "pval",
                "padj",
                "odds.ratio",
                "neg.log2.odds",
                "Z"
            )) %>%
            write.csv(paste0(study, "--overall--arm-gain-mutexclusivity.csv"),
                row.names = FALSE,
                na = "",
                quote = TRUE
            )
    }
}
