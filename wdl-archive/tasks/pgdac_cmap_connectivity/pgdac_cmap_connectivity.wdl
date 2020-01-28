task pgdac_cmap_connectivity {
  File tarball
  Int permutations
  Int? rankpt_n
  Int? rankpt_threshold
  Array[File] subset_scores
  Array[File]? permutation_scores
  String scores_dir = "cmap-subset-scores"
  String permutation_dir = "cmap-permutation-scores"
  String outFile = "pgdac_cmap-output.tar"

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions

  command {
    set -euo pipefail
    # setup output scores directory ...
    if [ ! -d ${scores_dir} ]; then 
      mkdir ${scores_dir} 
    fi
    # ... and copy subset scores
    mv ${sep=" " subset_scores} ${scores_dir}
    
    # same for perumations scores ...
    if [ ${permutations} -gt 0 ]; then
      if [ ! -d ${permutation_dir} ]; then
        mkdir ${permutation_dir}
      fi
      mv ${sep=" " permutation_scores} ${permutation_dir}
    fi
    
    # combine shards/gather and run conectivity score calculations
    /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh CMAPconn -i ${tarball} -o ${outFile} -CMAPscr ${scores_dir} -CMAPnperm ${permutations} -CMAPpmt ${permutation_dir} -CMAPrpt ${default="4" rankpt_n} -CMAPth ${default="85" rankpt_threshold}
  }

  output {
    File outputs = "${outFile}"
  }

  runtime {
    docker : "broadcptac/pgdac_main:2"
    memory : select_first ([memory, 32]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 64]) + " SSD"
    cpu : select_first ([num_threads, permutations+1]) + ""
    preemptible : select_first ([num_preemptions, 0])
  }

  meta {
    author : "D. R. Mani"
    email : "manidr@broadinstitute.org"
  }
}

workflow pgdac_cmap_connectivity_workflow {
	call pgdac_cmap_connectivity
}
