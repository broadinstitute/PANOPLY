task pgdac_association {
  File tarball   # output from pgdac_harmonize or pgdac_normalize_ms_data
  File? groupsFile
  String type
  String? subType
  File? params
  String codeDir = "/prot/proteomics/Projects/PGDAC/src"
  String outFile = "pgdac_association-output.tar"

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions


  command {
    set -euo pipefail
    /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh assoc -i ${tarball} -t ${type} -c ${codeDir} -o ${outFile} ${"-g " + groupsFile} ${"-m " + subType} ${"-p " + params}
  }

  output {
    File outputs = "${outFile}"
  }

  runtime {
    docker : "broadcptac/pgdac_main:1"
    memory : select_first ([memory, 4]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 5]) + " SSD"
    cpu : select_first ([num_threads, 1]) + ""
    preemptible : select_first ([num_preemtions, 0])
  }

  meta {
    author : "D. R. Mani"
    email : "manidr@broadinstitute.org"
  }
}

task pgdac_cluster {
  File tarball   # output from pgdac_harmonize or pgdac_normalize_ms_data
  String type
  String? subType
  File? params
  String codeDir = "/prot/proteomics/Projects/PGDAC/src"
  String outFile = "pgdac_cluster-output.tar"

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions


  command {
    set -euo pipefail
    /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh cluster -i ${tarball} -t ${type} -c ${codeDir} -o ${outFile} ${"-m " + subType} ${"-p " + params}
  }

  output {
    File outputs = "${outFile}"
  }

  runtime {
    docker : "broadcptac/pgdac_main:1"
    memory : select_first ([memory, 4]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 5]) + " SSD"
    cpu : select_first ([num_threads, 1]) + ""
    preemptible : select_first ([num_preemtions, 0])
  }

  meta {
    author : "D. R. Mani"
    email : "manidr@broadinstitute.org"
  }
}

task pgdac_cna_correlation {
  File tarball   # output from pgdac_cna_setup
  String type
  String? subType
  File? params
  String outFile = "pgdac_cna_correlation-output.tar"

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions


  command {
    set -euo pipefail
    /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh CNAcorr -i ${tarball} -t ${type} -o ${outFile} ${"-m " + subType} ${"-p " + params}
  }

  output {
    File outputs = "${outFile}"
  }

  runtime {
    docker : "broadcptac/pgdac_main:1"
    memory : select_first ([memory, 4]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 5]) + " SSD"
    cpu : select_first ([num_threads, 1]) + ""
    preemptible : select_first ([num_preemtions, 0])
  }

  meta {
    author : "D. R. Mani"
    email : "manidr@broadinstitute.org"
  }
}

task pgdac_cna_setup {
  File tarball   # output from pgdac_harmonize
  File? groupsFile
  String type
  String? subType
  File? params
  String codeDir = "/prot/proteomics/Projects/PGDAC/src"
  String outFile = "pgdac_cna_setup-output.tar"

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions


  command {
    set -euo pipefail
    /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh CNAsetup -i ${tarball} -t ${type} -c ${codeDir} -o ${outFile} ${"-g " + groupsFile} ${"-m " + subType} ${"-p " + params}
  }

  output {
    File outputs = "${outFile}"
  }

  runtime {
    docker : "broadcptac/pgdac_main:1"
    memory : select_first ([memory, 4]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 5]) + " SSD"
    cpu : select_first ([num_threads, 1]) + ""
    preemptible : select_first ([num_preemtions, 0])
  }

  meta {
    author : "D. R. Mani"
    email : "manidr@broadinstitute.org"
  }
}

task pgdac_harmonize {
  File tarball
  File rnaExpr
  File cnaExpr
  String type
  String? subType
  File? params
  String codeDir = "/prot/proteomics/Projects/PGDAC/src"
  String dataDir = "/prot/proteomics/Projects/PGDAC/data"
  String outFile = "pgdac_harmonize-output.tar"

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions


  command {
    set -euo pipefail
    /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh harmonize -i ${tarball} -t ${type} -c ${codeDir} -d ${dataDir} -rna ${rnaExpr} -cna ${cnaExpr} -o ${outFile} ${"-m " + subType} ${"-p " + params}
  }

  output {
    File outputs = "${outFile}"
  }

  runtime {
    docker : "broadcptac/pgdac_main:1"
    memory : select_first ([memory, 4]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 5]) + " SSD"
    cpu : select_first ([num_threads, 1]) + ""
    preemptible : select_first ([num_preemtions, 0])
  }

  meta {
    author : "D. R. Mani"
    email : "manidr@broadinstitute.org"
  }
}

task pgdac_normalize_ms_data {
  File tarball
  String type
  String? subType
  File? params
  String codeDir = "/prot/proteomics/Projects/PGDAC/src"
  String outFile = "pgdac_normalize_ms_data-output.tar"

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions

  command {
    set -euo pipefail
    /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh noramlize -i ${tarball} -t ${type} -c ${codeDir} -o ${outFile} ${"-m " + subType} ${"-p " + params}
  }

  output {
    File outputs = "${outFile}"
  }

  runtime {
    docker : "broadcptac/pgdac_main:1"
    memory : select_first ([memory, 4]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 5]) + " SSD"
    cpu : select_first ([num_threads, 1]) + ""
    preemptible : select_first ([num_preemtions, 0])
  }

  meta {
    author : "D. R. Mani"
    email : "manidr@broadinstitute.org"
  }
}

task pgdac_parse_sm_table {
  File SMtable
  File exptDesign
  String analysisDir
  String type
  String? subType
  File? params
  String codeDir = "/prot/proteomics/Projects/PGDAC/src"
  String dataDir = "/prot/proteomics/Projects/PGDAC/data"
  String outFile = "pgdac_parse_sm_table-output.tar"

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions

  command {
    set -euo pipefail
    /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh inputSM -s ${SMtable} -t ${type} -r ${analysisDir} -c ${codeDir} -d ${dataDir} -e ${exptDesign} -o ${outFile} ${"-m " + subType} ${"-p " + params}
  }

  output {
    File outputs = "${outFile}"
  }

  runtime {
    docker : "broadcptac/pgdac_main:2"
    memory : select_first ([memory, 4]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 5]) + " SSD"
    cpu : select_first ([num_threads, 1]) + ""
    preemptible : select_first ([num_preemtions, 0])
  }

  meta {
    author : "D. R. Mani"
    email : "manidr@broadinstitute.org"
  }
}

task pgdac_parse_sm_table_report {
  File tarball
  String label
  String type = "proteome"
  String tmpDir = "tmp"

  Int? memory
  Int? disk_space
  Int? num_threads

  command {
    set -euo pipefail
    Rscript /src/rmd-normalize.r ${tarball} ${label} ${type} ${tmpDir}
  }

  output {
    File report = "norm.html"
  }

  runtime {
    docker : "broadcptac/pgdac_cpdb:4"
    memory : select_first ([memory, 4]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 5]) + " SSD"
    cpu : select_first ([num_threads, 1]) + ""
  }

  meta {
    author : "Karsten Krug"
    email : "karsten@broadinstitute.org"
  }
}

task pgdac_rna_protein_correlation {
  File tarball
  File rnaExpr
  String type
  String? subType
  File? params
  String codeDir = "/prot/proteomics/Projects/PGDAC/src"
  String outFile = "pgdac_rna_protein_correlation-output.tar"

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions


  command {
    set -euo pipefail
    /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh RNAcorr -i ${tarball} -t ${type} -c ${codeDir} -rna ${rnaExpr} -o ${outFile} ${"-m " + subType} ${"-p " + params}
  }

  output {
    File outputs = "${outFile}"
  }

  runtime {
    docker : "broadcptac/pgdac_main:1"
    memory : select_first ([memory, 4]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 5]) + " SSD"
    cpu : select_first ([num_threads, 1]) + ""
    preemptible : select_first ([num_preemtions, 0])
  }

  meta {
    author : "D. R. Mani"
    email : "manidr@broadinstitute.org"
  }
}

task pgdac_rna_protein_correlation_report {
  File tarball
  String label
  String type = "proteome"
  String tmpDir = "tmp"
  Float fdr = 0.05

  Int? memory
  Int? disk_space
  Int? num_threads

  command {
    set -euo pipefail
    Rscript /src/rmd-rna-seq-correlation.r ${tarball} ${label} ${type} ${fdr} ${tmpDir}
  }

  output {
    File report = "rna-corr.html"
  }

  runtime {
    docker : "broadcptac/pgdac_cpdb:4"
    memory : select_first ([memory, 4]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 5]) + " SSD"
    cpu : select_first ([num_threads, 1]) + ""
  }

  meta {
    author : "Karsten Krug"
    email : "karsten@broadinstitute.org"
  }
}

task pgdac_sampleqc {
  File tarball   # output from pgdac_harmonize
  String type
  String? subType
  File? params
  String codeDir = "/prot/proteomics/Projects/PGDAC/src"
  String outFile = "pgdac_sampleqc-output.tar"

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions


  command {
    set -euo pipefail
    /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh sampleQC -i ${tarball} -t ${type} -c ${codeDir} -o ${outFile} ${"-m " + subType} ${"-p " + params}
  }

  output {
    File outputs = "${outFile}"
  }

  runtime {
    docker : "broadcptac/pgdac_main:1"
    memory : select_first ([memory, 4]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 5]) + " SSD"
    cpu : select_first ([num_threads, 1]) + ""
    preemptible : select_first ([num_preemtions, 0])
  }

  meta {
    author : "D. R. Mani"
    email : "manidr@broadinstitute.org"
  }
}

workflow pgdac_main_pipeline {
  File SMtable
  File exptDesign
  File rnaData
  File cnaData
  String analysisDir
  Float corr_fdr
  String? dataType="proteome"
  String? dataSubType
  File? additionalParameters
  File? cna_groups
  File? association_groups

  call pgdac_parse_sm_table {
    input:
      SMtable=SMtable, 
      analysisDir=analysisDir, 
      exptDesign=exptDesign,
      type=dataType,
      subType=dataSubType,
      params=additionalParameters
  }

  call pgdac_normalize_ms_data {
    input:
      tarball=pgdac_parse_sm_table.outputs,
      type=dataType,
      subType=dataSubType,
      params=additionalParameters
  }

  call pgdac_parse_sm_table_report {
    input: 
      tarball=pgdac_normalize_ms_data.outputs,
      label=analysisDir
  }

  call pgdac_rna_protein_correlation {
    input: 
      tarball=pgdac_normalize_ms_data.outputs, 
      rnaExpr=rnaData,
      type=dataType,
      subType=dataSubType,
      params=additionalParameters
  }
  
  call pgdac_rna_protein_correlation_report { 
    input:
      tarball=pgdac_rna_protein_correlation.outputs,
      label=analysisDir,
      fdr=corr_fdr
  }

  call pgdac_harmonize {
    input: 
      tarball=pgdac_rna_protein_correlation.outputs, 
      rnaExpr=rnaData,
      cnaExpr=cnaData,
      type=dataType,
      subType=dataSubType,
      params=additionalParameters
 }

  call pgdac_sampleqc {
    input:
      tarball=pgdac_harmonize.outputs,
      type=dataType,
      subType=dataSubType,
      params=additionalParameters
  }

  call pgdac_cna_setup {
    input: 
      tarball=pgdac_sampleqc.outputs, 
      groupsFile=cna_groups,
      type=dataType,
      subType=dataSubType,
      params=additionalParameters
  }

  call pgdac_cna_correlation {
    input:
      tarball=pgdac_cna_setup.outputs,
      type=dataType,
      subType=dataSubType,
      params=additionalParameters
  }

  call pgdac_association {
    input: 
      tarball=pgdac_cna_correlation.outputs, 
      groupsFile=association_groups,
      type=dataType,
      subType=dataSubType,
      params=additionalParameters
  }

  call pgdac_cluster {
    input:
      tarball=pgdac_association.outputs,
      type=dataType,
      subType=dataSubType,
      params=additionalParameters
  }

  output {
    File output=pgdac_cluster.outputs
    File norm_report=pgdac_parse_sm_table_report.report
    File corr_report=pgdac_rna_protein_correlation_report.report
  }}

