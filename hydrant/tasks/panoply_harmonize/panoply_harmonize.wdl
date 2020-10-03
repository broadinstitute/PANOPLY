#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
task panoply_harmonize {
  File inputData
  File rnaExpr
  File cnaExpr
  String type
  String standalone
  String? analysisDir
  String? subType
  File yaml
  String? pomeGeneIdCol
  String? cnaGeneIdCol
  String? rnaGeneIdCol
  Float? na_max
  String? duplicate_gene_policy
  String? gene_id_col


  String codeDir = "/prot/proteomics/Projects/PGDAC/src"
  String dataDir = "/prot/proteomics/Projects/PGDAC/data"
  String outFile = "panoply_harmonize-output.tar"

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions


  command {
    set -euo pipefail
    Rscript /prot/proteomics/Projects/PGDAC/src/parameter_manager.r \
    --module harmonize \
    --master_yaml ${yaml} \
    ${"--pome_gene_id_col " + pomeGeneIdCol} \
    ${"--cna_gene_id_col " + cnaGeneIdCol} \
    ${"--rna_gene_id_col " + rnaGeneIdCol} \
    ${"--na_max " + na_max} \
    ${"--duplicate_gene_policy " + duplicate_gene_policy} \
    ${"--gene_id_col " + gene_id_col}

    if [[ ${standalone} = false ]]; then
      /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh harmonize \
                  -i ${inputData} \
                  -t ${type} \
                  -c ${codeDir} \
                  -d ${dataDir} \
                  -rna ${rnaExpr} \
                  -cna ${cnaExpr} \
                  -o ${outFile} \
                  ${"-m " + subType} \
                  -p "config-custom.r";
    else
      /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh harmonize \
                  -f ${inputData} \
                  -r ${analysisDir} \
                  -t ${type} \
                  -c ${codeDir} \
                  -d ${dataDir} \
                  -rna ${rnaExpr} \
                  -cna ${cnaExpr} \
                  -o ${outFile} \
                  ${"-m " + subType} \
                  -p "config-custom.r";
    fi
  }

  output {
    File outputs = "${outFile}"
  }

  runtime {
    docker : "broadcptacdev/panoply_harmonize:latest"
    memory : select_first ([memory, 12]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 20]) + " SSD"
    cpu : select_first ([num_threads, 1]) + ""
    preemptible : select_first ([num_preemptions, 0])
  }

  meta {
    author : "Ramani Kothadia"
    email : "proteogenomics@broadinstitute.org"
  }
}

workflow panoply_harmonize_workflow {
    String standalone
    File inputData
    File rnaExpr
    File cnaExpr
    String dataType
    String? analysisDir
    File yaml
    Float? na_max
    String? duplicate_gene_policy
    String? gene_id_col

  call panoply_harmonize {
    input:
      inputData=inputData,
      rnaExpr=rnaExpr,
      cnaExpr=cnaExpr,
      analysisDir=analysisDir,
      standalone=standalone,
      type=dataType,
      yaml=yaml,
      na_max=na_max,
      duplicate_gene_policy=duplicate_gene_policy,
      gene_id_col=gene_id_col
  }
}
