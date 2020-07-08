task panoply_unified_pre
{
  File yaml_config

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions

  command {
    set -euo pipefail
    Rscript /prot/proteomics/Projects/PGDAC/src/unified_pre.r -c ${yaml_config}
  }

  output {
    File file_of_types = "types.txt"
  } 

  runtime {
    docker : "broadcptacdev/panoply_unified_pre:latest"
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

workflow panoply_unified_pre_workflow {
  call panoply_unified_pre
  Array[File] array_of_types = read_lines(panoply_unified_pre.file_of_types)
  output {
    Array[File] array_output = array_of_types
  }
}
