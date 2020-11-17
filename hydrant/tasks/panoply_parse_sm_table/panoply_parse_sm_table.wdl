#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
task panoply_parse_sm_table {
  File SMtable
  File exptDesign
  String analysisDir
  String type
  File yaml
  String? labelType
  String? speciesFilter
  Int? ndigits
  String outFile = "panoply_parse_sm_table-output.tar"

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions

  command {
    set -euo pipefail
    codeDir="/prot/proteomics/Projects/PGDAC/src"
    dataDir="/prot/proteomics/Projects/PGDAC/data"
    Rscript /prot/proteomics/Projects/PGDAC/src/parameter_manager.r \
    --module parse_sm_table \
    --master_yaml ${yaml} \
    ${"--label_type " + labelType} \
    ${"--species_filter " + speciesFilter} \
    ${"--ndigits " + ndigits}
    /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh inputSM -s ${SMtable} -t ${type} -r ${analysisDir} -c $codeDir -d $dataDir -e ${exptDesign} -o ${outFile} -p "/prot/proteomics/Projects/PGDAC/src/new-config-custom.r" -y "final_output_params.yaml"
  }

  output {
    File outputs = "${outFile}"
  }

  runtime {
    docker : "broadcptacdev/panoply_parse_sm_table:latest"
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


workflow panoply_parse_sm_table_workflow {
	call panoply_parse_sm_table
}
