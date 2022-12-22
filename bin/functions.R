p.to.Z <- function(p) {
    #' Converts a two-tailed p-value to standardized Z score
    if (!is.numeric(p)) {
        p <- as.numeric(p)
    }
    if (!is.numeric(p)) {
        stop("p was not coercible to numeric type")
    }
    ll <- length(which(p < 0 | p > 1))
    if (ll > 0) {
        warning(ll, " invalid p-values set to NA")
        p[p < 0 | p > 1] <- NA
    }
    O <- qnorm((p / 2), F)
    O[!is.finite(O)] <- NA
    return(-O)
}
coxph.zscores <- function(fit) {
    #' Given a `survival` fit object, extract the Z-scores of variables
    with(fit, coefficients / sqrt(diag(var)))
}
mutmatrix <- function(neither, gene.not.cna,
                      cna.not.gene, both) {
    # `` Produces a 2x2 matrix in format:
    # ``       cna  !cna
    # `` gene  ...  ...
    # ``!gene  ...  ...
    matrix(
        data = c(
            both, gene.not.cna,
            cna.not.gene, neither
        ),
        ncol = 2,
        byrow = TRUE,
        dimnames = list(c("gene", "!gene"), c("cna", "!cna"))
    )
}
mutexclusivity <- function(geneMatrix, cnaMatrix, minimumPatientsWcna = 100) {
    #' Given a binary matrix of mutated genes x samples and
    #' binary matrix of cna x samples, calculate the exlusivity of
    #' mutated gene and cna across samples
    #'
    #' On a per gene x cna comparions, NA samples are dropped.

    # Assert that the matrices are in the proper type
    if (typeof(geneMatrix) != typeof(TRUE)) {
        stop("geneMatrix is not a logical matrix")
    }
    if (typeof(cnaMatrix) != typeof(TRUE)) {
        stop("cnaMatrix is not a logical matrix")
    }

    # Ensure we compare each sample to itself across mutations and CNAs
    .samples <- sort(intersect(rownames(geneMatrix), rownames(cnaMatrix)))
    geneMatrix <- geneMatrix[.samples, ]
    cnaMatrix <- cnaMatrix[.samples, ]

    # Keep only CNAs with enough patients
    cnaMatrix <- cnaMatrix[, colSums(cnaMatrix, na.rm = TRUE) >= minimumPatientsWcna]

    genes <- colnames(geneMatrix)
    cnas <- colnames(cnaMatrix)
    .matrix <- matrix(
        nrow = length(genes),
        ncol = length(cnas),
        dimnames = list(genes, cnas)
    )
    matrices <- list(
        neither = .matrix,
        both = .matrix,
        gene.not.cna = .matrix,
        cna.not.gene = .matrix,
        fisher.p = .matrix,
        fisher.odds = .matrix
    )

    # Note: Order of loops has been checked for cache efficiency using `microbenchmark::microbenchmark(...)`
    for (gene in genes) {
        for (cna in cnas) {
            g <- geneMatrix[, gene]
            c <- cnaMatrix[, cna]

            # Remove any samples which are NA in either vector
            g.notNA <- rownames(geneMatrix)[which(!is.na(g))]
            c.notNA <- rownames(cnaMatrix)[which(!is.na(c))]
            .samples <- sort(intersect(g.notNA, c.notNA))
            g <- g[.samples]
            c <- c[.samples]

            # `table` is not used because it does not always produce a 2x2 matrix
            # if either c or g contain only one level it is 2x1 or 1x2
            both <- sum(g & c)
            gene.not.cna <- sum(g & !c)
            cna.not.gene <- sum(!g & c)
            neither <- sum(!(g | c))

            matrices$neither[gene, cna] <- neither
            matrices$both[gene, cna] <- both
            matrices$gene.not.cna[gene, cna] <- gene.not.cna
            matrices$cna.not.gene[gene, cna] <- cna.not.gene

            fisher <- fisher.test(mutmatrix(
                neither = neither,
                gene.not.cna = gene.not.cna,
                cna.not.gene = cna.not.gene,
                both = both
            ))
            matrices$fisher.p[gene, cna] <- fisher$p.value
            matrices$fisher.odds[gene, cna] <- fisher$estimate
        }
    }

    matrices
}

find_skip <- function(file, pattern, n = 5) {
    min(grep(pattern, readr::read_lines(file = file, n_max = n))) - 1
}
