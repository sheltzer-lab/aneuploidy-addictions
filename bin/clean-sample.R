#!/usr/bin/env Rscript

library(assertr, quietly = TRUE)
library(tidyverse, quietly = TRUE)

source(here::here("bin/functions.R"))

Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 5)

args <- commandArgs(trailingOnly = TRUE)

samples <- lapply(unlist(stringr::str_split(args[1], ",")), function(fh) {
    read_tsv(fh,
        skip = find_skip(fh, "PATIENT_ID"),
        col_types = cols(.default = "c")
    ) %>%
        select(PATIENT_ID, SAMPLE_ID, CANCER_TYPE, SAMPLE_TYPE)
}) %>% bind_rows()

# Determine which samples to keep as representative of a patient
# We only want to keep one sample per patient, but some patients have
# limited (CNA) data
samplesWCnaData <- lapply(unlist(stringr::str_split(args[2], ",")), function(fh) {
    read_tsv(fh,
        skip = find_skip(fh, "Hugo_Symbol"),
        col_types = cols(.default = "c")
    ) %>%
        select(-any_of(c("Hugo_Symbol", "Entrez_Gene_Id", "Cytoband"))) %>%
        colnames()
}) %>%
    unlist() %>%
    unique()
patientsWMultipleCancers <- samples %>%
    group_by(PATIENT_ID) %>%
    count(CANCER_TYPE) %>%
    filter(n > 1) %>%
    pluck("PATIENT_ID")
preferredSample <- samples %>%
    select(PATIENT_ID, SAMPLE_ID, SAMPLE_TYPE) %>%
    # Only consider samples with CNA data
    filter(SAMPLE_ID %in% samplesWCnaData) %>%
    # Remove patients with multiple cancers
    filter(!(PATIENT_ID %in% patientsWMultipleCancers)) %>%
    # Order by SAMPLE_TYPE so we can prioritize
    # Preference: ... < Metastasis < Recurrence < Primary
    mutate(SAMPLE_TYPE = factor(SAMPLE_TYPE, ordered = TRUE) %>%
        fct_relevel("Metastasis", "Recurrence", "Primary", after = Inf)) %>%
    group_by(PATIENT_ID) %>%
    arrange(desc(SAMPLE_TYPE)) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    pluck("SAMPLE_ID")

samples %>%
    select(PATIENT_ID, SAMPLE_ID, CANCER_TYPE) %>%
    filter(SAMPLE_ID %in% preferredSample) %>%
    verify(has_all_names("PATIENT_ID", "SAMPLE_ID", "CANCER_TYPE")) %>%
    verify(SAMPLE_ID %in% preferredSample) %>%
    assert(is_uniq, PATIENT_ID, SAMPLE_ID) %>%
    write.csv("clean_sample.csv", row.names = FALSE, na = "", quote = TRUE)
