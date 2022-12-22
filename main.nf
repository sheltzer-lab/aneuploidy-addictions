#!/usr/bin/env nextflow
nextflow.enable.dsl=2 

// The cBioPortal study (IDs) to analyze (comma-separated string)
params.studies = "msk_met_2021"

// Minimum number of patients in a cancer type in order to continue investigating that cancer (integer)
params.patientsCutoff = 25

// The percentage of patients with mutatation in a gene in order to continue investigating that gene (decimal)
params.percentMutated = 0.02

// Minimum number of patients with arm gain (integer)
params.patientsCutoffCna = 100

// Percent of patients with arm gain, used when patientsCutoffCna is not exceeded (decimal)
params.patientsCutoffCnaRelative = 0.02

// The limited gene set to investigate; adjusts for hypermutated cancers (CSV file with hugoGeneSymbol column)
params.geneSet = "${projectDir}/ref/impact-505.csv"

// The length of chromosomal arms (CSV file)
params.armDefinition = "${projectDir}/ref/hg19-chromosome-arm-lengths.csv"

workflow {
  def study = Channel.fromList(params.studies?.tokenize(','))
  def patientsCutoff = Channel.value(params.patientsCutoff)
  def percentMutated = Channel.value(params.percentMutated)
  def patientsCutoffCna = Channel.value(params.patientsCutoffCna)
  def patientsCutoffCnaRelative = Channel.value(params.patientsCutoffCnaRelative)
  def geneSet = Channel.fromPath(params.geneSet, checkIfExists: true) \
                | splitCsv(header: true, quote: '"') \
                | map { row -> row."hugoGeneSymbol" } \
                | collect
                | map { it -> it.join(",") }
  def armDefinition = Channel.fromPath(params.armDefinition, checkIfExists: true).first()

  // Download HGNC
  HGNC(armDefinition, geneSet)

  // Download study
  GIT_CLONE_DATAHUB()
  GIT_LFS_STUDY(GIT_CLONE_DATAHUB.out, study)

  // Clean data
  CLEAN_SAMPLE(GIT_LFS_STUDY.out.sample.groupTuple() | join(GIT_LFS_STUDY.out.cna.groupTuple()))
  CLEAN_PATIENT(GIT_LFS_STUDY.out.patient.groupTuple() | join(CLEAN_SAMPLE.out.clean), patientsCutoff)
  CLEAN_MUTATION(GIT_LFS_STUDY.out.mutation.groupTuple() | join(CLEAN_PATIENT.out.clean), geneSet, percentMutated)
  CLEAN_SCORE(GIT_LFS_STUDY.out.score.groupTuple() | join(CLEAN_PATIENT.out.clean))

  // Process data
  MUTUAL_EXCLUSIVITY(CLEAN_MUTATION.out.clean \
                     | join(CLEAN_SCORE.out.clean) \
                     | join(CLEAN_PATIENT.out.clean), 
                     patientsCutoffCna,
                     patientsCutoffCnaRelative)

  // Produce figures and tables for publication
  VIZ_FIGURES(
    MUTUAL_EXCLUSIVITY.out.map{it -> [it[0], "arm-gain", it[1].grep( ~/.*gain.*/ )]} \
    | join(CLEAN_MUTATION.out.clean ) \
    | join(CLEAN_SCORE.out.clean) \
    | join(CLEAN_PATIENT.out.clean)
  )
  TAB_SUPPLEMENTAL(
    MUTUAL_EXCLUSIVITY.out.map{it -> [it[0], "arm-gain", it[1].grep( ~/.*gain.*/ )]} \
    | join(CLEAN_MUTATION.out.clean ) \
    | join(CLEAN_SCORE.out.clean) \
    | join(CLEAN_PATIENT.out.clean)
  )
}

process HGNC {
  conda "${projectDir}/env/HGNC.yml"
  label 'process_medium'
  
  input:
    path definition
    val geneset

  output:
    path('hgnc.csv')

  script:
    """
    hgnc.R ${definition} '${geneset}'
    """
}

process GIT_CLONE_DATAHUB {
  conda "${projectDir}/env/git-lfs.yml"
  label 'process_short'
  publishDir "${projectDir}"

  output:
    path 'datahub'

  script:
    """
    git lfs install --skip-repo --skip-smudge
    git clone https://github.com/cBioPortal/datahub.git
    cd datahub
    git lfs install --local --skip-smudge
    """
}

process GIT_LFS_STUDY {
  maxForks 1  // Helps with rate limiting
  conda "${projectDir}/env/git-lfs.yml"
  label 'process_short'
  tag "${study}"

  input:
    path datahub
    val study

  output:
    tuple val(study), path("${study}_data_clinical_patient.txt"), emit: patient
    tuple val(study), path("${study}_data_clinical_sample.txt"), emit: sample
    tuple val(study), path("${study}_data_mutations.txt"), emit: mutation
    tuple val(study), path("${study}_data_cna.txt"), emit: cna
    tuple val(study), path("${study}_scores.csv"), emit: score

  script:
    """
    cd ${datahub}
    git lfs pull -I public/${study}/data_clinical_patient.txt
    git lfs pull -I public/${study}/data_clinical_sample.txt
    git lfs pull -I public/${study}/data_mutations.txt
    git lfs pull -I public/${study}/data_cna.txt
    
    cd -
    ln -s ${datahub}/public/${study}/data_clinical_patient.txt ${study}_data_clinical_patient.txt
    ln -s ${datahub}/public/${study}/data_clinical_sample.txt ${study}_data_clinical_sample.txt
    ln -s ${datahub}/public/${study}/data_mutations.txt ${study}_data_mutations.txt
    ln -s ${datahub}/public/${study}/data_cna.txt ${study}_data_cna.txt
    ln -s ${projectDir}/ref/${study}_scores.csv ${study}_scores.csv
    """
}

process CLEAN_PATIENT {
  conda "${projectDir}/env/cleaning.yml"
  label 'process_medium'
  tag "${study}"

  input:
    tuple val(study), 
          path(patientList), 
          path(sample)
    val patientsCutoff

  output:
    tuple val(study), path("clean_patient.csv"), emit: clean

  script:
    def patients =  patientList.join(",")
    """
    clean-patient.R '${patients}' ${sample} ${patientsCutoff}
    """
}

process CLEAN_SAMPLE {
  conda "${projectDir}/env/cleaning.yml"
  label 'process_medium'
  tag "${study}"

  input:
    tuple val(study), 
          path(sampleList),
          path(cnaList)

  output:
    tuple val(study), path("clean_sample.csv"), emit: clean

  script:
    def samples = sampleList.join(",")
    def cnas = cnaList.join(",")
    """
    clean-sample.R '${samples}' '${cnas}'
    """
}

process CLEAN_MUTATION {
  conda "${projectDir}/env/cleaning.yml"
  label 'process_medium'
  tag "${study}"

  input:
    tuple val(study), 
          path(mutationList), 
          path(patient)
    val geneset
    val percentMutated

  output:
    tuple val(study), path("clean_mutation.csv"), emit: clean

  script:
    def mutations =  mutationList.join(",")
    """
    clean-mutation.R '${mutations}' ${patient} '${geneset}' '${percentMutated}'
    """
}

process CLEAN_SCORE {
  conda "${projectDir}/env/cleaning.yml"
  label 'process_medium'
  tag "${study}"

  input:
    tuple val(study), 
          path(scoreList), 
          path(patient)

  output:
    tuple val(study), path("clean_score.csv"), emit: clean

  script:
    def scores =  scoreList.join(",")
    """
    clean-score.R '${scores}' ${patient}
    """
}

process MUTUAL_EXCLUSIVITY {
  conda "${projectDir}/env/mutual-exclusivity.yml"
  label 'process_long'
  tag "${study}"
  publishDir "results/${study}/mutexclusivity/", mode: 'copy', overwrite: true

  input:
    tuple val(study),
          path(mutation),
          path(score),
          path(patient)
    val patientsCutoffCna
    val patientsCutoffCnaRelative

  output:
    tuple val(study),
          path('*-mutexclusivity.csv')

  script:
    """
    mutual-exclusivity.R ${study} ${mutation} ${score} ${patientsCutoffCna} ${patient} ${patientsCutoffCnaRelative}
    """
}

process VIZ_FIGURES {
  conda "${projectDir}/env/figures.yml"
  label 'process_medium'
  tag "${study}-${type}"
  publishDir "results/${study}/viz/${type}/figures/", mode: 'copy', overwrite: true

  input:
    tuple val(study),
          val(type),
          path(mutexList),
          path(mutation),
          path(score),
          path(patient)

  output:
    path('*.png')

  script:
    def mutexs = mutexList.join(",")
    """
    viz-figures.R ${study} ${type} ${mutation} ${score} ${patient} ${mutexs}
    """
}

process TAB_SUPPLEMENTAL {
  conda "${projectDir}/env/tables.yml"
  label 'process_medium'
  tag "${study}-${type}"
  publishDir "results/${study}/tab/${type}/supplemental/", mode: 'copy', overwrite: true

  input:
    tuple val(study),
          val(type),
          path(mutexList),
          path(mutation),
          path(score),
          path(patient)

  output:
    path('Supplemental-Z-Tables.xlsx')

  script:
    def mutexs = mutexList.join(",")
    """
    tab-supplemental.R ${study} ${type} ${mutation} ${score} ${patient} ${mutexs}
    """
}
