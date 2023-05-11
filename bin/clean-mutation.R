#!/usr/bin/env Rscript

library(assertr, quietly = TRUE)
library(tidyverse, quietly = TRUE)

# Determine the path to local shared R definitions
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)

source(file.path(script.basename, "functions.R"))
source(file.path(script.basename, "constants.R"))

args <- commandArgs(trailingOnly = TRUE)

patient <- read.csv(args[2], stringsAsFactors = FALSE) %>%
    verify(has_all_names("PATIENT_ID", "STATUS", "MONTHS", "SAMPLE_ID", "CANCER_TYPE")) %>%
    assert(is_uniq, PATIENT_ID, SAMPLE_ID) %>%
    select(SAMPLE_ID, CANCER_TYPE)
geneset <- unlist(stringr::str_split(args[3], ","))
percentMutated <- as.numeric(args[4])

# Read data and limit to:
#  1. The mutations of interest
#  2. The samples of interest
lapply(unlist(stringr::str_split(args[1], ",")), function(fh) {
    read_tsv(fh,
        skip = find_skip(fh, "Hugo_Symbol"),
        col_types = cols(.default = "c")
    ) %>%
        select(Hugo_Symbol, Tumor_Sample_Barcode, Variant_Classification) %>%
        # Only include mutations that affect the protein.
        filter(!(Variant_Classification %in% ignoreMutations))
}) %>%
    bind_rows() %>%
    # Keep only those samples with a match in cleaned samples
    inner_join(patient, by = c("Tumor_Sample_Barcode" = "SAMPLE_ID")) %>%
    # Determine the HUGO symbols which we will keep.
    # By percentage mutated within a cancer and by geneset
    # Calculate the mutational cutoff for each cancer type
    group_by(CANCER_TYPE) %>%
    mutate(NPatients = length(unique(Tumor_Sample_Barcode))) %>%
    mutate(CUTOFF = percentMutated * NPatients) %>%
    # Count the number of patients in a cancer type with a mutation in a gene
    group_by(CANCER_TYPE, Hugo_Symbol) %>%
    mutate(NMutated = length(unique(Tumor_Sample_Barcode))) %>%
    ungroup() %>%
    # Remove low-mutation symbols
    filter(CUTOFF <= NMutated) %>%
    # Remove symbols not in gene set
    filter(Hugo_Symbol %in% geneset) %>%
    # Prep file and write
    select(Hugo_Symbol, Tumor_Sample_Barcode, Variant_Classification) %>%
    verify(has_all_names("Hugo_Symbol", "Tumor_Sample_Barcode", "Variant_Classification")) %>%
    assert(in_set(geneset), Hugo_Symbol) %>%
    verify(Tumor_Sample_Barcode %in% patient$SAMPLE_ID) %>%
    assert(in_set(ignoreMutations, inverse = TRUE), Variant_Classification) %>%
    write.csv("clean_mutation.csv", row.names = FALSE, na = "", quote = TRUE)
