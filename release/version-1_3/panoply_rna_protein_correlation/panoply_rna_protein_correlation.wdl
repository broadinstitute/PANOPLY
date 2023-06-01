#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
task panoply_rna_protein_correlation {
  File inputData
  File rnaExpr
  String type
  String standalone
  File yaml
  Int? rnaSDthreshold
  Int? profilePlotTopN
  String? analysisDir
  String outFile = "panoply_rna_protein_correlation-output.tar"

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions


  command {
    set -euo pipefail
    codeDir="/prot/proteomics/Projects/PGDAC/src"
    dataDir="/prot/proteomics/Projects/PGDAC/data"

    Rscript /prot/proteomics/Projects/PGDAC/src/parameter_manager.r \
    --module rna_protein_correlation \
    --master_yaml ${yaml} \
    ${"--rna_sd_threshold " + rnaSDthreshold} \
    ${"--profile_plot_top_n " + profilePlotTopN}
    if [[ ${standalone} = false ]]; then
      /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh RNAcorr \
                  -i ${inputData} \
                  -t ${type} \
                  -c $codeDir \
                  -rna ${rnaExpr} \
                  -o ${outFile} \
                  -p "/prot/proteomics/Projects/PGDAC/src/new-config-custom.r" \
                  -y "final_output_params.yaml";
    else
      /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh RNAcorr \
                  -f ${inputData} \
                  -t ${type} \
                  -c $codeDir \
                  -d $dataDir \
                  -rna ${rnaExpr} \
                  -r ${analysisDir} \
                  -o ${outFile} \
                  -p "/prot/proteomics/Projects/PGDAC/src/new-config-custom.r" \
                  -y "final_output_params.yaml";
    fi
  }

  output {
    File outputs = "${outFile}"
  }

  runtime {
    docker : "broadcptac/panoply_rna_protein_correlation:1_3"
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

workflow panoply_rna_protein_correlation_workflow {
  File rnaExpr
  String dataType
  File inputData
  String standalone
  String? analysisDir
  Int? rnaSDthreshold
  Int? profilePlotTopN

  call panoply_rna_protein_correlation {
    input:
      inputData=inputData,
      type=dataType,
      rnaExpr=rnaExpr,
      analysisDir=analysisDir,
      standalone=standalone,
      rnaSDthreshold=rnaSDthreshold,
      profilePlotTopN=profilePlotTopN
  } 
}