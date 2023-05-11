[![DOI](https://zenodo.org/badge/581325839.svg)](https://zenodo.org/badge/latestdoi/581325839)

# Aneuploidy Addictions

A script which conducts the computational analysis associated with the Sheltzer Lab's investigation into aneuploidy addictions.

## Scores

Aneuploidy scores follow the method of: https://www.nature.com/articles/s41467-020-14286-0

## Dependencies

This workflow has been designed such that only `nextflow` and `conda` should be necessary -- all other dependencies are setup as needed.

```console
$ nextflow -version

    N E X T F L O W
    version 22.04.5 build 5708
    created 15-07-2022 16:09 UTC (12:09 EDT)
    cite doi:10.1038/nbt.3820
    http://nextflow.io

```

```console
$ conda --version
conda 4.14.0
```

It should be sufficient to run the following on a machine with both of the above programs. This will first download the workflow from GitHub (with associated reference data), then run the workflow building conda environments as needed.

```console
nextflow run sheltzer-lab/aneuploidy-addictions -with-conda
```

## Workflow

Rather than repeat what can be read in the Nextflow script itself, the description here is meant to provide a "plain language" non-computational overview of the workflow. I will intermix searchable code elements in `monospaced font` so that one may cross-reference the description here with the code directly.

### Purpose

The primary purpose of this analysis is to complement the laboratory work done in the Sheltzer Lab by analyzing public human cancer data available from the [cBioPortal](https://www.cbioportal.org/). In order to do this, we must select which cBioPortal `--studies` we wish to analyze, i.e. `msk_met_2021`. There is a requirement that this study also have arm gain data in `ref/`.

### Parameters

Beyond `--studies`, workflow parameters include: `--patientsCutoff`, how many patients must be in a cancer type in order to analyze it; `--percentMutated`, what percentage of patients must have a mutation in a gene in order to analyze mutation data of that gene; `--patientsCutoffCna`, how many patients must present with gain/loss in order to continue analyzing the relevant CNA data; `--patientsCutoffCnaRelative`, sometimes `--patientsCutoffCna` is too stringent and we must loosen our requirements for rarer cancer types so this parameter allows us to treat the CNA data similar to mutation data via a relative cutoff of what percentage of patients must present with CNA in order to analyze CNA data; `--geneSet`, limits the genes we investigate to those with known cancer implication, reducing the high number of false positives which are an emergent property of genome-wide analysis; lastly `--armDefinition`, which allows us to map the genomic position to the chromosomal arm based on a given reference genome.

### Cleaning

All data must be cleaned. The procedure used here cleans all incoming data to initially minimizing the amount of data flowing through the system while also ensuring consistent sample selection. Some samples do not have complete data so we select a representative sample for each patient to ensure one sample per patient -- this procedure is done in the `CLEAN_SAMPLE` process -- which is then used to filter all other incoming data. Representative samples are those that have CNA/genomic data, while also not being from a patient with multiple cancers. Where a patient has multiple CNA/genomic profiles, we prioritize selecting the primary sample, followed by recurrence, and lastly metastasis. Beyond sample selection, this is also the stage at which gene selection is performed. In order to control for hypermutation we limit the investigated genes to those in the IMPACT-505 geneset and remove any genes with low levels of mutation.

## MISC

This project also included SMASH data processing. The code for processing SMASH sequencing data can be found at [sheltzer-lab/smash](https://github.com/sheltzer-lab/smash). The code for building standardized SMASH results plots can be found at [sheltzer-lab/smash_viz](https://github.com/sheltzer-lab/smash_viz).
