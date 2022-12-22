#!/usr/bin/env Rscript

library(assertr, quietly = TRUE)
library(tidyverse, quietly = TRUE)

source(here::here("bin/functions.R"))
source(here::here("bin/constants.R"))

args <- commandArgs(trailingOnly = TRUE)

preferredSample <- read.csv(args[2], stringsAsFactors = FALSE) %>%
    verify(has_all_names("PATIENT_ID", "STATUS", "MONTHS", "SAMPLE_ID", "CANCER_TYPE")) %>%
    assert(is_uniq, PATIENT_ID, SAMPLE_ID) %>%
    pluck("SAMPLE_ID") %>%
    unique()

lapply(unlist(stringr::str_split(args[1], ",")), function(fh) {
    read_csv(fh, col_types = cols(.default = "c")) %>%
        mutate(across(!c(PATIENT_ID, SAMPLE_ID), as.numeric)) %>%
        filter(SAMPLE_ID %in% preferredSample)
}) %>%
    bind_rows() %>%
    distinct() %>%
    verify(has_all_names("PATIENT_ID", "SAMPLE_ID")) %>%
    verify(SAMPLE_ID %in% preferredSample) %>%
    assert(is_uniq, PATIENT_ID, SAMPLE_ID) %>%
    assert(within_bounds(-1, 1), any_of(chroms)) %>%
    write.csv("clean_score.csv", row.names = FALSE, na = "", quote = TRUE)
