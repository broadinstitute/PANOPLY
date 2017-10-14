

task cna_analysis {
  File rna
  File cna
  File pome
  String prefix
  Int jidMax
  Int jid
  String codeDir = "/prot/proteomics/Projects/PGDAC/src"

  Int? memory
  Int? disk_space
  Int? num_threads

  command {
    set -euo pipefail
    # setup directories and code
    cp ${codeDir}/cna-analysis.r ${codeDir}/generate-cna-plots.r .
    if [ ! -d ${prefix}-output ]; then 
      mkdir ${prefix}-output 
    fi
    # run cna analysis for corresponding shard / gather
    Rscript cna-analysis.r ${jid} ${jidMax} ${prefix} ${rna} ${cna} ${pome}
  }

  output {
    File rna_cna_corr = "${prefix}-output/mrna-vs-cna-corr${jid}.csv"
    File rna_cna_pval = "${prefix}-output/mrna-vs-cna-pval${jid}.csv"
    File pome_cna_corr = "${prefix}-output/pome-vs-cna-corr${jid}.csv"
    File pome_cna_pval = "${prefix}-output/pome-vs-cna-pval${jid}.csv"
  }

  runtime {
    docker : "broadcptac/pgdac_basic:1"
    memory : select_first ([memory, 4]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 8]) + " SSD"
    cpu : select_first ([num_threads, 1]) + ""

  }

  meta {
    author : "D. R. Mani"
    email : "manidr@broadinstitute.org"
  }
}



task gather_results_and_plot {
  String prefix
  Int jidMax
  Array[File] rna_vs_cna_corr
  Array[File] rna_vs_cna_pval
  Array[File] pome_vs_cna_corr
  Array[File] pome_vs_cna_pval
  String codeDir = "/prot/proteomics/Projects/PGDAC/src"
  String dataDir = "/prot/proteomics/Projects/PGDAC/data"
  
  Int? memory
  Int? disk_space
  Int? num_threads

  command {
    set -euo pipefail
    # setup directories and code
    cp ${codeDir}/cna-analysis.r ${codeDir}/generate-cna-plots.r .
    cp ${dataDir}/chr-length.csv ${dataDir}/gene-location.csv .
    if [ ! -d ${prefix}-output ]; then 
      mkdir ${prefix}-output 
    fi
    # copy results from scatter operation
    mv ${sep=" " rna_vs_cna_corr} ${prefix}-output
    mv ${sep=" " rna_vs_cna_pval} ${prefix}-output
    mv ${sep=" " pome_vs_cna_corr} ${prefix}-output
    mv ${sep=" " pome_vs_cna_pval} ${prefix}-output
    # run cna analysis for corresponding shard / gather
    Rscript cna-analysis.r 0 ${jidMax} ${prefix} NULL NULL NULL
  }

  output {
    Array[File] tables=glob ("${prefix}-*-vs-*.csv")
    File plot="${prefix}-cna-plot.png"
  }

  runtime {
    docker : "broadcptac/pgdac_basic:1"
    memory : select_first ([memory, 16]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 16]) + " SSD"
    cpu : select_first ([num_threads, 1]) + ""
  }

  meta {
    author : "D. R. Mani"
    email : "manidr@broadinstitute.org"
  }
}




workflow run_cna_analysis {
  File rna
  File cna
  File pome
  String prefix
  File jidsFile
  Array[Int] jids = read_lines ("${jidsFile}")
  Int jidMax = length (jids)
  

  scatter (i in jids) {
    call cna_analysis {
      input:
        rna=rna,
        cna=cna,
        pome=pome,
        prefix=prefix,
        jidMax=jidMax,
        jid=i
    }
  }
  
  call gather_results_and_plot {
    input:
      prefix=prefix,
      jidMax=jidMax,
      rna_vs_cna_corr=cna_analysis.rna_cna_corr,
      rna_vs_cna_pval=cna_analysis.rna_cna_pval,
      pome_vs_cna_corr=cna_analysis.pome_cna_corr,
      pome_vs_cna_pval=cna_analysis.pome_cna_pval
  }

  output {
    Array[File] tables = gather_results_and_plot.tables
    File plot = gather_results_and_plot.plot
  }
}

