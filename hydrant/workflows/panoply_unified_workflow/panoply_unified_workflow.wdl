#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_normalize_ms_data_workflow/versions/1/plain-WDL/descriptor" as normalize_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_main/versions/14/plain-WDL/descriptor" as main_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_blacksheep_workflow/versions/3/plain-WDL/descriptor" as blacksheep_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_so_nmf_workflow/versions/18/plain-WDL/descriptor" as so_nmf_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_mo_nmf_gct/versions/5/plain-WDL/descriptor" as mo_nmf_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_immune_analysis_workflow/versions/3/plain-WDL/descriptor" as immune_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_make_pairs_workflow/versions/3/plain-WDL/descriptor" as make_pairs_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_unified_assemble_results/versions/5/plain-WDL/descriptor" as assemble_wdl


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
  String run_nmf #'true' or 'false'

  # Normalize specific optional params:
  String? normalizeProteomics # "true" or "false"

  
  Array[File] omes = ["${prote_ome}", "${phospho_ome}", "${acetyl_ome}", "${ubiquityl_ome}"]
  Array[File] data = ["${rna_data}", "${cna_data}"]
  
  # Turn the type data (RNA+CNA) into an array of pairs:
  call make_pairs_wdl.panoply_make_pairs_workflow as type_pairs {
    input:
      files = data,
      suffix = "-aggregate.gct"

  }

  # Make the array of ome pairs to input into normalize:
  call make_pairs_wdl.panoply_make_pairs_workflow as ome_pairs {
    input:
      files = omes,
      suffix = "-aggregate.gct"

  }

  ### NORMALIZE:
  ### Normalize the data first so downstream modules (NMF etc) can run in parallel to main:
  scatter (pair in ome_pairs.zipped) {
    if ("${pair.right}" != "") {
      call normalize_wdl.panoply_normalize_ms_data_workflow as norm {
        input:
          input_pome="${pair.right}",
          ome_type="${pair.left}",
          job_identifier="${job_id}-${pair.left}",
          yaml="${yaml}",
          normalizeProteomics=normalizeProteomics
      }
    }
  }

  # Make an array of pairs from the normalized data for input to main and other modules:
  call make_pairs_wdl.panoply_make_pairs_workflow as norm_pairs {
    input:
      files = norm.normalized_data_table,
      suffix = "-normalized_table-output.gct"

  }

  ### MAIN:
  scatter (pair in norm_pairs.zipped) {
    if ("${pair.right}" != "") {
        call main_wdl.panoply_main as pome {
          input:
            ## include all required arguments from above
            input_pome="${pair.right}",
            ome_type="${pair.left}",
            job_identifier="${job_id}-${pair.left}",
            run_ptmsea="${run_ptmsea}",
            run_cmap = "${run_cmap}",
            run_nmf = "false",
            input_cna="${cna_data}",
            input_rna_v3="${rna_data}",
            sample_annotation="${sample_annotation}",
            yaml="${yaml}"
        }
        
      
    }
  }
  

  # This takes the array of pairs of normalized proteomics data and combines it with the array of pairs of RNA+CNA data for Blacksheep use:
  Array[Pair[String?,File?]] all_pairs = flatten([norm_pairs.zipped,type_pairs.zipped])
  
  ### BLACKSHEEP:
  scatter (pair in all_pairs) {
    if ("${pair.right}" != "") {
      call blacksheep_wdl.panoply_blacksheep_workflow as outlier {
        input:
          input_gct = "${pair.right}",
          master_yaml = "${yaml}",
          output_prefix = "${pair.left}",
          type = "${pair.left}"
      }
    }
  }
  
  ### Single-ome NMF
  call so_nmf_wdl.panoply_so_nmf_workflow as so_nmf {
    input:
      yaml = yaml,
      job_id = job_id,
      prote_ome = norm.normalized_data_table[0],
      phospho_ome = norm.normalized_data_table[1],
      acetyl_ome = norm.normalized_data_table[2],
      ubiquityl_ome = norm.normalized_data_table[3],
      rna_data = rna_data,
      cna_data = cna_data
  }

  ### Multi-omics NMF:
  if ( run_nmf == "true" ){
    call mo_nmf_wdl.panoply_mo_nmf_gct_workflow as mo_nmf {
      input:
        yaml_file = yaml,
        label = job_id,
        omes = norm.normalized_data_table,
        rna_ome = rna_data,
        cna_ome = cna_data
    }
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
      norm_report = norm.normalize_report,
      rna_corr_report = pome.rna_corr_report,
      cna_corr_report = pome.cna_corr_report,
      omicsev_report = pome.omicsev_report,
      sampleqc_report = pome.sample_qc_report,
      assoc_report = pome.association_report,
      blacksheep_tar = outlier.blacksheep_tar,
      blacksheep_report = outlier.blacksheep_report,
      so_nmf_results = so_nmf.nmf_results,
      so_nmf_reports = so_nmf.nmf_reports,
      mo_nmf_tar = mo_nmf.nmf_clust,
      mo_nmf_report = mo_nmf.nmf_clust_report,
      mo_nmf_ssgsea_tar = mo_nmf.nmf_ssgsea,
      mo_nmf_ssgsea_report = mo_nmf.nmf_ssgsea_report,
      immune_tar = immune.outputs,
      immune_report = immune.report

  }
  
  output {
    File all_results = panoply_unified_assemble_results.all_results
    File all_reports = panoply_unified_assemble_results.all_reports
  }
 }
