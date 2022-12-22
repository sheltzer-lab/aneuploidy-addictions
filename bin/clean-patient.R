#!/usr/bin/env Rscript

library(assertr, quietly = TRUE)
library(tidyverse, quietly = TRUE)
library(lubridate, quietly = TRUE)

source(here::here("bin/functions.R"))

args <- commandArgs(trailingOnly = TRUE)

patients <- lapply(unlist(stringr::str_split(args[1], ",")), function(fh) {
    read_tsv(fh,
        skip = find_skip(fh, "PATIENT_ID"),
        col_types = cols(.default = "c")
    )
}) %>%
    bind_rows()

samples <- read.csv(args[2], stringsAsFactors = FALSE) %>%
    verify(has_all_names("PATIENT_ID", "SAMPLE_ID", "CANCER_TYPE")) %>%
    assert(is_uniq, PATIENT_ID, SAMPLE_ID)

nPatientCutoff <- as.numeric(args[3])

# Determine which STATUS/MONTHS data exists
# Preference is in written order so RFS_STATUS is selected over OS_STATUS if the former exists
status.name <- detect(colnames(patients), ~ . %in% c("RFS_STATUS", "OS_STATUS"))
months.name <- detect(colnames(patients), ~ . %in% c("RFS_MONTHS", "OS_MONTHS"))

preferredSample <- unique(samples$SAMPLE_ID)

tmp <- patients %>%
    select(
        PATIENT_ID,
        STATUS = all_of(status.name),
        MONTHS = all_of(months.name),
        any_of("CANCER_TYPE_ACRONYM")
    ) %>%
    mutate(
        STATUS = str_detect(STATUS, "1"), # String containing 0/1 for status
        MONTHS = dmonths(as.numeric(MONTHS)) %/% dmonths(1)
    ) %>%
    left_join(samples %>% select(SAMPLE_ID, PATIENT_ID, CANCER_TYPE),
        by = c("PATIENT_ID" = "PATIENT_ID")
    ) %>%
    {
        if ("CANCER_TYPE_ACRONYM" %in% names(.)) select(., -CANCER_TYPE) %>% rename(CANCER_TYPE = CANCER_TYPE_ACRONYM) else .
    } %>%
    filter(SAMPLE_ID %in% preferredSample)

cancerTypes <- tmp %>%
    group_by(CANCER_TYPE) %>%
    summarise(N = length(unique(SAMPLE_ID))) %>%
    # Remove cancer types w/o enough patients
    filter(nPatientCutoff <= N) %>%
    pluck("CANCER_TYPE") %>%
    unique()

tmp %>%
    filter(CANCER_TYPE %in% cancerTypes) %>%
    verify(has_all_names("PATIENT_ID", "STATUS", "MONTHS", "SAMPLE_ID", "CANCER_TYPE")) %>%
    verify(SAMPLE_ID %in% preferredSample) %>%
    assert(is_uniq, PATIENT_ID, SAMPLE_ID) %>%
    assert(in_set(TRUE, FALSE), STATUS) %>%
    verify(MONTHS >= 0 | is.na(MONTHS)) %>%
    write.csv("clean_patient.csv", row.names = FALSE, na = "", quote = TRUE)
