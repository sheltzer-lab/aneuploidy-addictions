chroms <- expand.grid(c("p", "q"), c(seq(1, 22, 1), "X", "Y")) %>%
    apply(1, function(row) paste0(row[[2]], row[[1]]))

ignoreMutations <- c(
    "Silent",
    "Intron",
    "5'UTR",
    "3'Flank",
    "5'Flank",
    "3'UTR",
    "IGR",
    "RNA"
)

genesOfFocus <- c(
    "APC",
    "BRAF",
    "CDKN2A",
    "EGFR",
    "KRAS",
    "MYC",
    "PIK3CA",
    "PTEN",
    "TP53"
)
