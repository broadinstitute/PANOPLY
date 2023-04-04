#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
task panoply_association {
  File inputData
  String type
  String standalone
  String? analysisDir
  File groupsFile
  File yaml
  Float? fdr_assoc
  Float? sample_na_max
  Float? nmiss_factor
  String? duplicate_gene_policy
  String? gene_id_col

  String outFile = "panoply_association-output.tar"

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions

  command {
    set -euo pipefail
    codeDir="/prot/proteomics/Projects/PGDAC/src"
    Rscript /prot/proteomics/Projects/PGDAC/src/parameter_manager.r \
    --module association \
    --master_yaml ${yaml} \
    ${"--fdr_assoc " + fdr_assoc} \
    ${"--sample_na_max " + sample_na_max} \
    ${"--nmiss_factor " + nmiss_factor} \
    ${"--duplicate_gene_policy " + duplicate_gene_policy} \
    ${"--gene_id_col " + gene_id_col}
    if [[ ${standalone} = false ]]; then
      /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh assoc \
                  -i ${inputData} \
                  -t ${type} \
                  -c $codeDir \
                  -o ${outFile} \
                  -g ${groupsFile} \
                  -y "final_output_params.yaml" \
                  -p "/prot/proteomics/Projects/PGDAC/src/new-config-custom.r";
    else
      /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh assoc \
                  -f ${inputData} \
                  -t ${type} \
                  -c $codeDir \
                  -r ${analysisDir} \
                  -o ${outFile} \
                  -g ${groupsFile} \
                  -y "final_output_params.yaml" \
                  -p "/prot/proteomics/Projects/PGDAC/src/new-config-custom.r"
    fi
  }

  output {
    File outputs = "${outFile}"
  }

  runtime {
    docker : "broadcptacdev/panoply_association:latest"
    memory : select_first ([memory, 16]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 40]) + " SSD"
    cpu : select_first ([num_threads, 1]) + ""
    preemptible : select_first ([num_preemptions, 0])
  }

  meta {
    author : "Ramani Kothadia"
    email : "proteogenomics@broadinstitute.org"
  }
}

workflow panoply_association_workflow {
  String standalone
  File inputData
  String? analysisDir
  File? groupsFile
  String dataType
  File yaml
  Float? fdr_assoc
  Float? sample_na_max
  Float? nmiss_factor
  String? duplicate_gene_policy
  String? gene_id_col

    
  call panoply_association {
    input:
      inputData=inputData,
      type=dataType,
      standalone=standalone,
      analysisDir=analysisDir,
      groupsFile=groupsFile,
      yaml=yaml,
      fdr_assoc=fdr_assoc,
      sample_na_max=sample_na_max,
      nmiss_factor=nmiss_factor,
      duplicate_gene_policy=duplicate_gene_policy,
      gene_id_col=gene_id_col
  }
}