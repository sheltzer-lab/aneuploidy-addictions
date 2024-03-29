#!/usr/bin/env Rscript
library(reshape2)
library(tidyverse)
library(openxlsx)

# Determine the path to local shared R definitions
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)

source(file.path(script.basename, "constants.R")) # Load `chroms`

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

cancers <- sort(unique(cancerMutexs$cancer))
cancerDFs <- lapply(cancers, function(c) {
    cancerMutexs %>%
        filter(cancer == c) %>%
        acast(gene ~ arm, value.var = "Z", fill = NA) %>%
        as.data.frame() %>%
        rownames_to_column("Gene")
}) %>%
    setNames(cancers)

overallDF <- overallMutex %>%
    acast(gene ~ arm, value.var = "Z", fill = NA) %>%
    as.data.frame() %>%
    rownames_to_column("Gene")

allTabs <- c(cancerDFs, list(overall = overallDF))

names(allTabs) <- names(allTabs) %>% str_trunc(width = 30)

write.xlsx(x = allTabs, file = "Supplemental-Z-Tables.xlsx", keepNA = TRUE)
