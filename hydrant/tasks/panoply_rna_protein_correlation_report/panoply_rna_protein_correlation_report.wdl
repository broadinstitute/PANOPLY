task panoply_rna_protein_correlation_report {
  File tarball
  File config_yaml
  
  String label
  String type
  String tmpDir
  #Float fdr

  Int? memory
  Int? disk_space
  Int? num_threads

  command {
    set -euo pipefail
    cfg=${config_yaml}
    echo "library(yaml);yaml=read_yaml('$cfg');fdr=yaml[['panoply_rna_protein_correlation_report']]\$fdr;writeLines(as.character(fdr), con='fdr.txt')" > cmd.r
    Rscript cmd.r
    fdr=`cat fdr.txt`
    Rscript /home/pgdac/src/rmd-rna-seq-correlation.r ${tarball} ${label} ${type} $fdr ${tmpDir}
  }

  output {
    File report = "rna-corr_" + label + ".html"
  }

  runtime {
    docker : "broadcptacdev/panoply_rna_protein_correlation_report:latest"
    memory : select_first ([memory, 8]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 20]) + " SSD"
    cpu : select_first ([num_threads, 1]) + ""
  }

  meta {
    author : "Karsten Krug"
    email : "karsten@broadinstitute.org"
  }
}


workflow panoply_rna_protein_correlation_report_workflow {
	call panoply_rna_protein_correlation_report
}
