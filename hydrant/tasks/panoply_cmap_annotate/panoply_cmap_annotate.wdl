#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
task panoply_cmap_annotate {
  File tarball                  # output from pgdac_cmap_connectivity
  File cmap_data_file           # CMAP level 5 geneKD data (gctx)
  File? cmap_enrichment_groups   # groups file (ala experiment design file)
  File yaml
  String? cmap_grp
  String? cmap_typ
  String outFile = "panoply_cmap-annotate-output.tar"

  Float? cna_threshold
  String? log_transform
  String? alpha

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions

  String cmap_group = "${if defined (cmap_grp) then cmap_grp else 'all'}"
  String cmap_type = "${if defined (cmap_typ) then cmap_typ else 'pome'}"

  command {
    set -euo pipefail
    Rscript /prot/proteomics/Projects/PGDAC/src/parameter_manager.r \
    --module cmap_analysis \
    --master_yaml ${yaml} \
    ${"--cna_threshold " + cna_threshold} \
    ${"--log_transform " + log_transform} \
    ${"--alpha " + alpha}
    Rscript /prot/proteomics/Projects/PGDAC/src/cmap-annotate.R  ${tarball} ${cmap_data_file} ${cmap_group} ${cmap_type} ${cmap_enrichment_groups} ${outFile} "cmap-config-custom.r"
  }

  output {
    File outputs = "${outFile}"
    File gsea_input = "${cmap_group}-cmap-${cmap_type}-gsea-input.gct"
  }

  runtime {
    docker : "broadcptacdev/panoply_cmap_annotate:latest"
    memory : select_first ([memory, 32]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 64]) + " SSD"
    cpu : select_first ([num_threads, 1]) + ""
    preemptible : select_first ([num_preemptions, 0])
  }

  meta {
    author : "D R Mani"
    email : "proteogenomics@broadinstitute.org"
  }
}


workflow panoply_cmap_annotate_workflow {
  call panoply_cmap_annotate
}
