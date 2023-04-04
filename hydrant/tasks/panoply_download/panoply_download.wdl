#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
task panoply_download
{
  File association_tar
  File ssgsea_ome_tar
  File ssgsea_rna_tar
  File? ptmsea
  File? so_nmf_tar
  File? so_nmf_ssgsea_tar
  File? omicsev_tar
  File? cosmo_tar
  String output_prefix
  String analysisDir
  String summary_tar = "panoply_main_summary.tar"
  String full_tar = "panoply_main_full.tar"
  Array[File] ssgsea_assoc_tars
  String ssgsea_assoc_dir = "ssgsea_assoc"

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions

  command {
    set -euo pipefail

    if [ ! -d ${ssgsea_assoc_dir} ]; then
      mkdir ${ssgsea_assoc_dir}
    fi

    index=0
    for file in ${sep=' ' ssgsea_assoc_tars} ; do
      basefilename=$(basename $file)
      index=$((index+1))
      cp $file ${ssgsea_assoc_dir}/$basefilename-$index.tar;
    done


    /prot/proteomics/Projects/PGDAC/src/download.sh \
        -t ${association_tar} \
        -o ${ssgsea_ome_tar} \
        -r ${ssgsea_rna_tar} \
        -a ${analysisDir} \
        -s ${ssgsea_assoc_dir} \
        ${"-p" + ptmsea} \
        ${"-n" + so_nmf_tar} \
        ${"-m" + so_nmf_ssgsea_tar} \
        ${"-e" + omicsev_tar} \
        ${"-c" + cosmo_tar};
    mv ${summary_tar} ${output_prefix}-${summary_tar}
    mv ${full_tar} ${output_prefix}-${full_tar}
  }

  output {
    File summary = "${output_prefix}-${summary_tar}"
    File full = "${output_prefix}-${full_tar}"
  } 

  runtime {
    docker : "broadcptacdev/panoply_download:latest"
    memory : select_first ([memory, 12]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 20]) + " SSD"
    cpu : select_first ([num_threads, 1]) + ""
    preemptible : select_first ([num_preemptions, 0])
  }

  meta {
    author : "Ramani Kothadia"
    email : "proteogenomics@broadinstitute.org"
  }
}

workflow panoply_download_workflow {
	call panoply_download
}
