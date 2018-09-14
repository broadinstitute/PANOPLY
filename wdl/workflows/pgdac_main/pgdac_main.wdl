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
    docker : "broadcptac/pgdac_main:2"
    memory : select_first ([memory, 12]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 20]) + " SSD"
    cpu : select_first ([num_threads, 1]) + ""
    preemptible : select_first ([num_preemptions, 0])
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
    docker : "broadcptac/pgdac_main:2"
    memory : select_first ([memory, 12]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 20]) + " SSD"
    cpu : select_first ([num_threads, 1]) + ""
    preemptible : select_first ([num_preemptions, 0])
  }

  meta {
    author : "D. R. Mani"
    email : "manidr@broadinstitute.org"
  }
}

task pgdac_cna_correlation_report {
  File tarball
  String label
  String type
  String tmpDir
  Float fdr

  Int? memory
  Int? disk_space
  Int? num_threads

  command {
    set -euo pipefail
    Rscript /home/pgdac/src/rmd-cna-analysis.r ${tarball} ${label} ${type} ${fdr} ${tmpDir}
  }

  output {
    File report = "cna-analysis_" + label + ".html"
  }

  runtime {
    docker : "broadcptac/pgdac_rmd:3"
    memory : select_first ([memory, 8]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 20]) + " SSD"
    cpu : select_first ([num_threads, 1]) + ""
  }

  meta {
    author : "Karsten Krug"
    email : "karsten@broadinstitute.org"
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
    docker : "broadcptac/pgdac_main:2"
    memory : select_first ([memory, 12]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 20]) + " SSD"
    cpu : select_first ([num_threads, 1]) + ""
    preemptible : select_first ([num_preemptions, 0])
  }

  meta {
    author : "D. R. Mani"
    email : "manidr@broadinstitute.org"
  }
}

task pgdac_cons_clust {
  File tarball   # output from pgdac_harmonize or pgdac_normalize_ms_data
  String type
  File? groupsFile
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
    echo ${type}
    /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh cluster -i ${tarball} -t ${type} -c ${codeDir} -o ${outFile} ${"-m " + subType} ${"-p " + params} ${"-g " + groupsFile}
   
  }

  output {
    File outputs = "${outFile}"
  }

  runtime {
    docker : "broadcptac/pgdac_cons_clust:2"
    memory : select_first ([memory, 16]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 20]) + " SSD"
    cpu : select_first ([num_threads, 8]) + ""
    preemptible : select_first ([num_preemptions, 0])
  }

  meta {
    author : "D. R. Mani, Karsten Krug"
    email : "karsten@broadinstitute.org"
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
    docker : "broadcptac/pgdac_main:2"
    memory : select_first ([memory, 12]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 20]) + " SSD"
    cpu : select_first ([num_threads, 1]) + ""
    preemptible : select_first ([num_preemptions, 0])
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
    /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh normalize -i ${tarball} -t ${type} -c ${codeDir} -o ${outFile} ${"-m " + subType} ${"-p " + params}
  }

  output {
    File outputs = "${outFile}"
  }

  runtime {
    docker : "broadcptac/pgdac_main:2"
    memory : select_first ([memory, 12]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 20]) + " SSD"
    cpu : select_first ([num_threads, 1]) + ""
    preemptible : select_first ([num_preemptions, 0])
  }

  meta {
    author : "D. R. Mani"
    email : "manidr@broadinstitute.org"
  }
}

task pgdac_normalize_ms_data_report {
  File tarball
  String label
  String type
  String tmpDir

  Int? memory
  Int? disk_space
  Int? num_threads

  command {
    set -euo pipefail
    Rscript /home/pgdac/src/rmd-normalize.r ${tarball} ${label} ${type} ${tmpDir}
  }

  output {
    File report = "norm_" + label + ".html"
  }

  runtime {
    docker : "broadcptac/pgdac_rmd:3"
    memory : select_first ([memory, 8]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 20]) + " SSD"
    cpu : select_first ([num_threads, 1]) + ""
  }

  meta {
    author : "Karsten Krug"
    email : "karsten@broadinstitute.org"
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
    memory : select_first ([memory, 12]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 20]) + " SSD"
    cpu : select_first ([num_threads, 1]) + ""
    preemptible : select_first ([num_preemptions, 0])
  }

  meta {
    author : "D. R. Mani"
    email : "manidr@broadinstitute.org"
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
    docker : "broadcptac/pgdac_main:2"
    memory : select_first ([memory, 12]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 20]) + " SSD"
    cpu : select_first ([num_threads, 1]) + ""
    preemptible : select_first ([num_preemptions, 0])
  }

  meta {
    author : "D. R. Mani"
    email : "manidr@broadinstitute.org"
  }
}

task pgdac_rna_protein_correlation_report {
  File tarball
  String label
  String type
  String tmpDir
  Float fdr

  Int? memory
  Int? disk_space
  Int? num_threads

  command {
    set -euo pipefail
    Rscript /home/pgdac/src/rmd-rna-seq-correlation.r ${tarball} ${label} ${type} ${fdr} ${tmpDir}
  }

  output {
    File report = "rna-corr_" + label + ".html"
  }

  runtime {
    docker : "broadcptac/pgdac_rmd:3"
    memory : select_first ([memory, 8]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 20]) + " SSD"
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
    docker : "broadcptac/pgdac_main:2"
    memory : select_first ([memory, 12]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 20]) + " SSD"
    cpu : select_first ([num_threads, 1]) + ""
    preemptible : select_first ([num_preemptions, 0])
  }

  meta {
    author : "D. R. Mani"
    email : "manidr@broadinstitute.org"
  }
}

task pgdac_sampleqc_report {
  File tarball
  String label
  String type
  String tmpDir

  Int? memory
  Int? disk_space
  Int? num_threads

  command {
    set -euo pipefail
    Rscript /home/pgdac/src/rmd-sample-qc.r ${tarball} ${label} ${type} ${tmpDir}
  }

  output {
    File report = "sample-qc_" + label + ".html"
  }

  runtime {
    docker : "broadcptac/pgdac_rmd:3"
    memory : select_first ([memory, 8]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 20]) + " SSD"
    cpu : select_first ([num_threads, 1]) + ""
  }

  meta {
    author : "Karsten Krug"
    email : "karsten@broadinstitute.org"
  }
}

workflow pgdac_main_pipeline {
  File SMtable
  File exptDesign
  File rnaData
  File cnaData
  String analysisDir
  Float corr_fdr
  String dataType
  String? dataSubType
  File? additionalParameters
  File? cna_groups
  File? association_groups
  File? cluster_enrichment_groups

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

  call pgdac_normalize_ms_data_report {
    input: 
      tarball=pgdac_normalize_ms_data.outputs,
      label=analysisDir,
      type=dataType,
			tmpDir="tmp"
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
      fdr=corr_fdr,
      type=dataType,
			tmpDir="tmp"
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

  call pgdac_sampleqc_report {
    input:
			tarball=pgdac_sampleqc.outputs,
			type=dataType,
			label=analysisDir,
			tmpDir="tmp"
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

  call pgdac_cna_correlation_report {
    input:
			tarball=pgdac_cna_correlation.outputs,
			type=dataType,
			label=analysisDir,
			fdr=corr_fdr,
			tmpDir="tmp"
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
      groupsFile=cluster_enrichment_groups,
      subType=dataSubType,
      params=additionalParameters
  }

  output {
    File output=pgdac_cluster.outputs
    File norm_report=pgdac_normalize_ms_data_report.report
    File rna_corr_report=pgdac_rna_protein_correlation_report.report
		File cna_corr_report=pgdac_cna_correlation_report.report
		File sample_qc=pgdac_sampleqc_report.report
  }}

