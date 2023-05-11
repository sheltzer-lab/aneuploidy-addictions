#!/usr/bin/env Rscript

# biomaRt is used explictly below as it mangles the global env
library(assertr, quietly = TRUE)
library(tidyverse, quietly = TRUE)

# Determine the path to local shared R definitions
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)

source(file.path(script.basename, "constants.R")) # Load `chroms`

args <- commandArgs(trailingOnly = TRUE)

definition <- args[1]
geneset <- unlist(stringr::str_split(args[2], ","))

chrom.arms <- read_csv(definition, col_types = cols(.default = "c")) %>%
    verify(has_all_names("chrom", "arm", "length", "percent")) %>%
    mutate(length = as.numeric(length)) %>%
    mutate(chromarm = factor(paste0(chrom, arm),
        ordered = TRUE,
        levels = chroms
    )) %>%
    arrange(chromarm) %>%
    mutate(
        start = as.numeric(cumsum(lag(length, default = 0))),
        end = as.numeric(cumsum(length))
    ) %>%
    mutate(chrom = ifelse(chrom == "X", 23, chrom)) %>%
    mutate(chrom = ifelse(chrom == "Y", 24, chrom)) %>%
    mutate(chrom = as.numeric(chrom)) %>%
    mutate(col = ifelse((chrom %% 2) == 0, TRUE, FALSE)) %>%
    mutate(chrom = ifelse(chrom == 23, "X", chrom)) %>%
    mutate(chrom = ifelse(chrom == 24, "Y", chrom)) %>%
    mutate(chrom = as.character(chrom))

# Used explictly rather than loaded globally as it mangles the global env
biomaRt::getBM(
    attributes = c(
        "hgnc_symbol",
        "chromosome_name",
        "band",
        "start_position",
        "end_position"
    ),
    filters = list("with_hgnc" = TRUE),
    mart = biomaRt::useEnsembl(
        biomart = "ensembl",
        dataset = "hsapiens_gene_ensembl"
    )
) %>%
    filter(chromosome_name %in% as.character(c(1:22, "X", "Y"))) %>%
    mutate(across(c("start_position", "end_position"), as.numeric)) %>%
    assert(not_na, start_position, end_position) %>%
    mutate(mid_position = ((end_position - start_position) / 2) + start_position) %>%
    mutate(arm = ifelse(str_detect(band, "p"), "p", "q")) %>%
    left_join(
        chrom.arms %>%
            group_by(chrom) %>%
            mutate(chrstart = min(start)) %>%
            ungroup() %>%
            select(chrom, arm, chrstart),
        by = c("chromosome_name" = "chrom", "arm" = "arm")
    ) %>%
    assert(not_na, chrstart, mid_position) %>%
    mutate(absmid = chrstart + mid_position) %>%
    filter(hgnc_symbol %in% geneset) %>%
    assert(is_uniq, hgnc_symbol) %>%
    assert(in_set(geneset), hgnc_symbol) %>%
    verify(start_position < mid_position & mid_position < end_position) %>%
    assert(in_set("p", "q"), arm) %>%
    write.csv("hgnc.csv",
        row.names = FALSE,
        na = "",
        quote = TRUE
    )
