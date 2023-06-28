#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#

################################################
task panoply_so_nmf_assemble_drivers {

  File nmf_tar
  String ome_type

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions
  
  command {
    Rscript /home/pgdac/src/concatenate_nmf_gct.R ${nmf_tar} ${ome_type}
  }
  
  output {
    Pair[String, File] driver_pair=(ome_type, glob("*.gct")[0]) #pull ome_type and the concatenated gct
  }
  
  runtime {
    docker : "broadcptacdev/panoply_so_nmf_assemble_drivers:latest"
    memory : select_first ([memory, 16]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 32]) + " SSD"
    cpu : select_first ([num_threads, 1]) + ""
    preemptible : select_first ([num_preemptions, 0])
  }

  meta {
    author : "C. Williams"
    email : "proteogenomics@broadinstitute.org"
  }
}
