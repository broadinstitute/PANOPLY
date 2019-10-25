task pgdac_download
{
  File cons_clust_tar
  File ssgsea_ome_tar
  File ssgsea_rna_tar
  File? ptmsea
  String analysisDir
  String summary_tar = "pgdac_main_summary.tar"
  String full_tar = "pgdac_main_full.tar"
  Array[File] ssgsea_assoc_tars
  Array[File] ssgsea_clust_tars
  String ssgsea_assoc_dir = "ssgsea_assoc"
  String ssgsea_clust_dir = "ssgsea_clust"

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions

  command {
    set -euo pipefail

    if [ ! -d ${ssgsea_assoc_dir} ]; then
      mkdir ${ssgsea_assoc_dir}
    fi

    if [ ! -d ${ssgsea_clust_dir} ]; then
      mkdir ${ssgsea_clust_dir}
    fi

    index=0
    for file in ${sep=' ' ssgsea_assoc_tars} ; do
      basefilename=$(basename $file)
      index=$((index+1))
      cp $file ${ssgsea_assoc_dir}/$basefilename-$index.tar;
    done

    index=0
    for file in ${sep=' ' ssgsea_clust_tars} ; do
      basefilename=$(basename $file)
      index=$((index+1))
      cp $file ${ssgsea_clust_dir}/$basefilename-$index.tar;
    done

    /prot/proteomics/Projects/PGDAC/src/download.sh \
        -c ${cons_clust_tar} \
        -o ${ssgsea_ome_tar} \
        -r ${ssgsea_rna_tar} \
        -a ${analysisDir} \
        -s ${ssgsea_assoc_dir} \
        -l ${ssgsea_clust_dir} \
        ${"-p" + ptmsea};
  }

  output {
    File summary = "${summary_tar}"
    File full = "${full_tar}"
  } 

  runtime {
    docker : "broadcptac/pgdac_download:1"
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

workflow pgdac_download_workflow {
	call pgdac_download
}
