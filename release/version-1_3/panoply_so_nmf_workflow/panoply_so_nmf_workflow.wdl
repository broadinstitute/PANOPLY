#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_so_nmf_gct/versions/5/plain-WDL/descriptor" as so_nmf_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_so_nmf_assemble_results/versions/3/plain-WDL/descriptor" as assemble_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_so_nmf_sankey_workflow/versions/5/plain-WDL/descriptor" as so_nmf_sankey_wdl


workflow panoply_so_nmf_workflow {
  File gene_set_database
  File? prote_ome
  File? phospho_ome
  File? acetyl_ome
  File? ubiquityl_ome
  File? nglyco_ome
  File? rna_data      #version 1.3 only!
  File? cna_data
  File? multiomic_nmf_tar
  File yaml
  String job_id

  String run_sankey = "true"

  if (defined(prote_ome)) {
    Pair[String, File?] prote_pair = ("prot", prote_ome)
  }
  if (defined(phospho_ome)) {
    Pair[String, File?] phospho_pair = ("pSTY", phospho_ome)
  }
  if (defined(acetyl_ome)) {
    Pair[String, File?] acetyl_pair = ("acK", acetyl_ome)
  }
  if (defined(ubiquityl_ome)) {
    Pair[String, File?] ubiquityl_pair = ("ubK", ubiquityl_ome)
  }
  if (defined(nglyco_ome)) {
    Pair[String, File?] nglyco_pair = ("nglyco", nglyco_ome)
  }
  
  if (defined(rna_data)) {
    Pair[String, File?] rna_pair = ("RNA", rna_data)
  }
  if (defined(cna_data)) {
    Pair[String, File?] cna_pair = ("CNA", cna_data)
  }
  
  Array[Pair[String, File?]] ome_pairs = select_all([prote_pair, phospho_pair, acetyl_pair, ubiquityl_pair, nglyco_pair, rna_pair, cna_pair]) #select non-null inputs
  


  ### MAIN:
  scatter (pair in ome_pairs) {
    call so_nmf_wdl.panoply_so_nmf_gct_workflow as so_nmf {
      input:
      yaml_file = yaml,
      label = "${job_id}-${pair.left}",
      ome = "${pair.right}",
      ome_type = "${pair.left}",
      gene_set_database = gene_set_database
    }
  }
  
  ## assemble final output combining results from all NMFs
  call assemble_wdl.panoply_so_nmf_assemble_results as nmf_assemble {
    input:
      so_nmf_tar = so_nmf.nmf_clust,
      so_nmf_report = so_nmf.nmf_clust_report,
      so_nmf_ssgsea_tar = so_nmf.nmf_ssgsea,
      so_nmf_ssgsea_report = so_nmf.nmf_ssgsea_report
  }
  

  ## generate sankey diagrams and report ()
  if ( run_sankey == "true" ){
    call so_nmf_sankey_wdl.panoply_so_nmf_sankey_workflow as nmf_sankey {
      input:
        so_nmf_tar = nmf_assemble.nmf_results,
        mo_nmf_tar = multiomic_nmf_tar,
        label = "${job_id}",
    }
  }


  output {
    File nmf_results = nmf_assemble.nmf_results
    File nmf_reports = nmf_assemble.nmf_reports
    File? sankey_tar = nmf_sankey.sankey_tar
    File? sankey_report = nmf_sankey.sankey_report
  }
 }