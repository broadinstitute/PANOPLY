#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_main/versions/3/plain-WDL/descriptor" as main_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac_MM:panoply_blacksheep_workflow_MM/versions/3/plain-WDL/descriptor" as blacksheep_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_mo_nmf_gct/versions/16/plain-WDL/descriptor" as mo_nmf_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac_MM:panoply_immune_analysis_workflow_MM/versions/8/plain-WDL/descriptor" as immune_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac_MM:panoply_make_pairs_workflow_MM/versions/1/plain-WDL/descriptor" as make_pairs_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac_MM:panoply_unified_assemble_results_MM/versions/22/plain-WDL/descriptor" as assemble_wdl


workflow panoply_unified_workflow {
  File? prote_ome
  File? phospho_ome
  File? acetyl_ome
  File? ubiquityl_ome
  File? rna_data      #version 1.3 only!
  File? cna_data
  File? sample_annotation
  File yaml
  String job_id
  String run_ptmsea
  String run_cmap
  
  Array[File] omes = ["${prote_ome}", "${phospho_ome}", "${acetyl_ome}", "${ubiquityl_ome}"]
  Array[File] data = ["${rna_data}", "${cna_data}"]
  
  call make_pairs_wdl.panoply_make_pairs_workflow as type_pairs {
    input:
      files = data,
      suffix = "-aggregate.gct"

  }

  call make_pairs_wdl.panoply_make_pairs_workflow as ome_pairs {
    input:
      files = omes,
      suffix = "-aggregate.gct"

  }

  ### MAIN:
  scatter (pair in ome_pairs.zipped) {
    if ("${pair.right}" != "") {
        call main_wdl.panoply_main as pome {
          input:
            ## include all required arguments from above
            input_pome="${pair.right}",
            ome_type="${pair.left}",
            job_identifier="${job_id}-${pair.left}",
            run_ptmsea="${run_ptmsea}",
            run_cmap = "${run_cmap}",
            input_cna="${cna_data}",
            input_rna_v3="${rna_data}",
            sample_annotation="${sample_annotation}",
            yaml="${yaml}"
        }
        
      
    }
  }
  
  call make_pairs_wdl.panoply_make_pairs_workflow as norm_pairs {
    input:
      files = pome.normalized_data_table,
      suffix = "-normalized_table-output.gct"

  }

  Array[Pair[String?,File?]] all_pairs = flatten([norm_pairs.zipped,type_pairs.zipped])
  
  ### BLACKSHEEP:
  scatter (pair in all_pairs) {
    if ("${pair.right}" != "") {
      if ("${pair.left}" != "cna") {
        call blacksheep_wdl.panoply_blacksheep_workflow as outlier {
          input:
            input_gct = "${pair.right}",
            master_yaml = "${yaml}",
            output_prefix = "${pair.left}",
            type = "${pair.left}"
        }
      }
    }
  }
  
  ### NMF:
  call mo_nmf_wdl.panoply_mo_nmf_gct_workflow as nmf {
    input:
      yaml_file = yaml,
      label = job_id,
      omes = pome.normalized_data_table,
      rna_ome = rna_data,
      cna_ome = cna_data
  }
  
  ### IMMUNE:
  if ( "${rna_data}" != '' ) {
    call immune_wdl.panoply_immune_analysis_workflow as immune {
      input:
          inputData=rna_data,
          standalone="true",
          type="rna",
          yaml=yaml,
          analysisDir=job_id,
          label=job_id
    }
  }
  
  ## assemble final output combining results from panoply_main, blacksheep immune_analysis and mo_nmf
  call assemble_wdl.panoply_unified_assemble_results {
    input:
      main_full = pome.panoply_full,
      main_summary = pome.summary_and_ssgsea,
      cmap_output = pome.cmap_output,
      cmap_ssgsea_output = pome.cmap_ssgsea_output,
      norm_report = pome.norm_report,
      rna_corr_report = pome.rna_corr_report,
      cna_corr_report = pome.cna_corr_report,
      sampleqc_report = pome.sample_qc_report,
      assoc_report = pome.association_report,
      cons_clust_report = pome.cons_clust_report,
      blacksheep_tar = outlier.blacksheep_tar,
      blacksheep_report = outlier.blacksheep_report,
      mo_nmf_tar = nmf.nmf_clust,
      mo_nmf_report = nmf.nmf_clust_report,
      mo_nmf_ssgsea_tar = nmf.nmf_ssgsea,
      mo_nmf_ssgsea_report = nmf.nmf_ssgsea_report,
      immune_tar = immune.outputs,
      immune_report = immune.report

  }
  
  output {
    File all_results = panoply_unified_assemble_results.all_results
    File all_reports = panoply_unified_assemble_results.all_reports
  }
 }