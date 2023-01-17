#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
task panoply_filter {
  File inputData
  String type
  String standalone
  String analysisDir
  File yaml
  String? filterProteomics
  Int? ndigits
  Float? naMax
  String? geneIdCol
  Float? sdFilterThreshold
  Float? minNumratioFraction
  Int? minNumratioProteome
  Int? minNumratioPTMs
  String? applySMfilter

  String outTar = "panoply_filter-output.tar"
  String outTable = "filter_table-output.gct"

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions

  command  <<<
    set -euo pipefail
    
    codeDir="/prot/proteomics/Projects/PGDAC/src"
    dataDir="/prot/proteomics/Projects/PGDAC/data"
    
    Rscript /prot/proteomics/Projects/PGDAC/src/parameter_manager.r \
      --module filter \
      --master_yaml ${yaml} \
      ${"--filter_proteomics " + filterProteomics} \
      ${"--ndigits " + ndigits} \
      ${"--gene_id_col " + geneIdCol} \
      ${"--na_max " + naMax} \
      ${"--sd_filter_threshold " + sdFilterThreshold} \
      ${"--min_numratio_fraction " + minNumratioFraction} \
      ${"--min_numratio_proteome " + minNumratioProteome} \
      ${"--min_numratio_ptms " + minNumratioPTMs} \
      ${"--apply_sm_filter " + applySMfilter}

    # Find the flag for filter.proteomics in the yaml:
    cfg='final_output_params.yaml'
    echo "library(yaml);yaml=read_yaml('$cfg');filt=yaml[['filter.proteomics']];writeLines(as.character(filt), con='filt.txt')" > cmd.r
    Rscript cmd.r
    filt=`cat filt.txt`
    
    # If not filtering input is the filtered tab and becomes the output, else run filtering:
    if [ $filt = FALSE ]; then
      cp ${inputData} ${type}-${outTable}
      tar -c -f ${outTar} ${type}-${outTable}
    else
      if [[ ${standalone} = false ]]; then
        /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh filter \
              -i ${inputData} \
              -t ${type} \
              -c $codeDir \
              -o ${outTar} \
              -p "/prot/proteomics/Projects/PGDAC/src/new-config-custom.r" \
              -y "final_output_params.yaml"
      else
        /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh filter \
              -a ${inputData} \
              -r ${analysisDir} \
              -t ${type} \
              -c $codeDir \
              -d $dataDir \
              -o ${outTar} \
              -p "/prot/proteomics/Projects/PGDAC/src/new-config-custom.r" \
              -y "final_output_params.yaml"
      fi
      # Grab the filtered gct to set as output with appropriate name
      outGCT=`find ${analysisDir}/filtered-data -type f -iname "*-ratio-filt-NArm.gct"`
      outTableName=${type}-${outTable} 
      cp $outGCT $outTableName
    fi
  >>>

  output {
    File outputs = "${type}-${outTable}"
    File output_tar = "${outTar}"
    File output_yaml = "final_output_params.yaml"
  }

  runtime {
    docker : "broadcptacdev/panoply_filter:latest"
    memory : select_first ([memory, 12]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 20]) + " SSD"
    cpu : select_first ([num_threads, 1]) + ""
    preemptible : select_first ([num_preemptions, 0])
  }

  meta {
    author : "C.M. Williams"
    email : "proteogenomics@broadinstitute.org"
  }
}

workflow panoply_filter {
  File inputData
  String dataType
  String standalone
  String analysisDir
  File yaml
  Int? ndigits
  String? geneIdCol
  Float? naMax
  Float? sdFilterThreshold
  Float? minNumratioFraction

  call panoply_filter {
    input:
      inputData=inputData,
      type=dataType,
      standalone=standalone,
      analysisDir=analysisDir,
      yaml=yaml,
      ndigits=ndigits,
      geneIdCol=geneIdCol,
      naMax=naMax,
      sdFilterThreshold=sdFilterThreshold,
      minNumratioFraction=minNumratioFraction
  }
}