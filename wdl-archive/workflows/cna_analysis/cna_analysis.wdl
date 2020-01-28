

task harmonize_data {
  File rna
  File cna
  File pome
  String analysisDir
  String codeDir = "/prot/proteomics/Projects/PGDAC/src"
  String dataDir = "/prot/proteomics/Projects/PGDAC/data"

  Int? memory
  Int? disk_space
  Int? num_threads

  command {
    set -euo pipefail
    # create matrix files from input gct -- harmonize both rows and columns
    /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh harmonize -c ${codeDir} -d ${dataDir} -r ${analysisDir} -rna ${rna} -cna ${cna} -f ${pome}
  }

  output {
    File outputs = "harmonize-output.tar"
  }

  runtime {
    docker : "broadcptac/pgdac_basic:1"
    memory : select_first ([memory, 8]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 16]) + " SSD"
    cpu : select_first ([num_threads, 1]) + ""
  }

  meta {
    author : "D. R. Mani"
    email : "manidr@broadinstitute.org"
  }
}



task cna_analysis_setup {
  File tarball
  File? groups    # expt-design-like file to use for subgroups
  String codeDir = "/prot/proteomics/Projects/PGDAC/src"

  Int? memory
  Int? disk_space
  Int? num_threads

  command {
    set -euo pipefail
    # setup directories and code (not needed to run, but for final tar file output); 
    # use matrix files from harmonization and create any subsets
    # create table of matrix files (tsv, one line per group, in order: rna, cna, pome)
    # determine subgroup list (separately, since this is a string and cannot be part of matrix files)
    # determine actual jidMax and create file with list of jid's
    /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh CNAsetup -i ${tarball} -c ${codeDir} ${"-g " + groups}
  }

  output {
    Array[Array[String]] matrixFiles = read_tsv ("file_table.tsv")
    Array[String] subgroups = read_lines ("subgroups.txt")
    File outputs = "CNAsetup-output.tar"
  }

  runtime {
    docker : "broadcptac/pgdac_basic:1"
    memory : select_first ([memory, 8]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 16]) + " SSD"
    cpu : select_first ([num_threads, 1]) + ""
  }

  meta {
    author : "D. R. Mani"
    email : "manidr@broadinstitute.org"
  }
}



task cna_analysis {
  File tarball
  String rna
  String cna
  String pome
  String analysisDir
  String prefix
  Int jidMax = 1
  Int jid = 1
  String codeDir = "/prot/proteomics/Projects/PGDAC/src"
  String dataDir = "/prot/proteomics/Projects/PGDAC/data"

  Int? memory
  Int? disk_space
  Int? num_threads

  command {
    set -euo pipefail
    # setup directories, code and supporting data files
    cp ${codeDir}/cna-analysis.r ${codeDir}/generate-cna-plots.r .
    cp ${dataDir}/gene-location.csv ${dataDir}/chr-length.csv .
    if [ ! -d ${prefix}-output ]; then 
      mkdir ${prefix}-output 
    fi
    # extract appropriate input files from tarball and copy to cwd
    tar -x -f ${tarball} ${analysisDir}/cna
    cp ${analysisDir}/cna/${rna} ${analysisDir}/cna/${cna} ${analysisDir}/cna/${pome} .
    # run cna analysis using a single processor (jid=jidMax=1)
    # with the WGCNA library this is very fast, and does not need a scatter
    Rscript cna-analysis.r ${jid} ${jidMax} ${prefix} ${rna} ${cna} ${pome}
    # collect results and plot
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



task assemble_results {
  Array[Array[File]] table_files
  Array[File] plot_files
  File tarball

  command {
    set -eu   # do not use -o pipefail -- results in error 141 when using ... | head -1 | ...
    # extract tarball in current directory and set $analysis_dir
    tar -x -f ${tarball}
    analysis_dir=`tar -t -f ${tarball} | head -1 | sed -e 's/\/.*//'`
    cd $analysis_dir
    # copy result tables/plots to appropriate location
    # (first flatten the table_files 2D array;
    #  using sep=" " creates ["item1.1", ... "item1.n"] ["item2.1", ... "item2.n"] ...
    #  flatten by removing [ ] , "
    file_list=`echo '${sep=" " table_files}' | tr -d '][,"'`
    cp $file_list cna
    cp ${sep=" " plot_files} cna
    # recreate new tarball for output
    cd ..
    tar -c -f CNA-output.tar $analysis_dir
  }

  output {
    File outputs = "CNA-output.tar"
  }

  runtime {
    docker : "broadcptac/pgdac_basic:1"
    memory : "4GB"
    disks : "local-disk 16 SSD"
    cpu : "1"
  }

  meta {
    author : "D. R. Mani"
    email : "manidr@broadinstitute.org"
  }
}




workflow run_cna_analysis_on_subgroups {
  File rna
  File cna
  File pome
  String analysisDir

  call harmonize_data {
    input:
      rna=rna,
      cna=cna,
      pome=pome,
      analysisDir=analysisDir
  }
  
  call cna_analysis_setup {
    input:
       tarball=harmonize_data.outputs
  }
  
  scatter (idx in range (length (cna_analysis_setup.subgroups))) {
    call cna_analysis as cna_s {
      input:
        tarball=cna_analysis_setup.outputs,
        analysisDir=analysisDir,
        prefix=cna_analysis_setup.subgroups[idx],
        rna=cna_analysis_setup.matrixFiles[idx][0],
        cna=cna_analysis_setup.matrixFiles[idx][1],
        pome=cna_analysis_setup.matrixFiles[idx][2]
    }
  }
  
  call assemble_results {
    input:
      table_files=cna_s.tables,
      plot_files=cna_s.plot,
      tarball=cna_analysis_setup.outputs
  }
  
  output {
    File final_output = assemble_results.outputs
  }
}

