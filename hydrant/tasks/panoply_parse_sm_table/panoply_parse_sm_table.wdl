#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
task panoply_parse_sm_table {
  File SMtable
  File exptDesign
  String analysisDir
  String type
  File yaml
  String? subType
  String? labelType
  String? applySMfilter
  String? speciesFilter
  Int? ndigits
  Float? naMax
  Float? sampleNaMax
  Float? minNumratioFraction
  Float? nmissFactor
  Float? sdFilterThreshold
  String? duplicateGenePolicy
  String codeDir = "/prot/proteomics/Projects/PGDAC/src"
  String dataDir = "/prot/proteomics/Projects/PGDAC/data"
  String outFile = "panoply_parse_sm_table-output.tar"

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions

  command {
    set -euo pipefail
    Rscript /prot/proteomics/Projects/PGDAC/src/parameter_manager.r \
    --module parse_sm_table \
    --master_yaml ${yaml} \
    ${"--label_type " + labelType} \
    ${"--apply_sm_filter " + applySMfilter} \
    ${"--species_filter " + speciesFilter} \
    ${"--ndigits " + ndigits} \
    ${"--na_max " + naMax} \
    ${"--sample_na_max " + sampleNaMax} \
    ${"--min_numratio_fraction " + minNumratioFraction} \
    ${"--nmiss_factor " + nmissFactor} \
    ${"--sd_filter_threshold " + sdFilterThreshold} \
    ${"--duplicate_gene_policy " + duplicateGenePolicy}
    /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh inputSM -s ${SMtable} -t ${type} -r ${analysisDir} -c ${codeDir} -d ${dataDir} -e ${exptDesign} -o ${outFile} ${"-m " + subType} -p "config-custom.r"
  }

  output {
    File outputs = "${outFile}"
  }

  runtime {
    docker : "broadcptac/panoply_parse_sm_table:dev"
    memory : select_first ([memory, 12]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 20]) + " SSD"
    cpu : select_first ([num_threads, 1]) + ""
    preemptible : select_first ([num_preemptions, 0])
  }

  meta {
    author : "Ramani Kothadia"
    email : "rkothadi@broadinstitute.org"
  }
}


workflow panoply_parse_sm_table_workflow {
	call panoply_parse_sm_table
}
