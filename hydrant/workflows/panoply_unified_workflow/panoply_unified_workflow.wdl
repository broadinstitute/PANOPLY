#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_select_all_pairs/versions/1/plain-WDL/descriptor" as select_pairs
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_normalize_filter_workflow/versions/14/plain-WDL/descriptor" as norm_filt_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_main/versions/52/plain-WDL/descriptor" as main_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_blacksheep_workflow/versions/13/plain-WDL/descriptor" as blacksheep_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_nmf_workflow/versions/37/plain-WDL/descriptor" as nmf_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_immune_analysis_workflow/versions/14/plain-WDL/descriptor" as immune_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_unified_assemble_results/versions/25/plain-WDL/descriptor" as assemble_wdl



workflow panoply_unified_workflow {
  File? prote_ome
  File? phospho_ome
  File? acetyl_ome
  File? ubiquityl_ome
  File? rna_data      #version 1.3 only!
  File? cna_data
  File yaml
  String job_id
  String run_cmap
  Boolean run_mo_nmf #'true' or 'false'
  Boolean run_so_nmf #'true' or 'false'
  String? run_ptmsea

  File groups_file
  File? groups_file_blacksheep
  File? groups_file_immune
  File? groups_file_nmf

  # Normalize specific optional params:
  String? normalizeProteomics # "true" or "false"
  String? filterProteomics # "true" or "false"

  # Organize omics data into pairs

  Array[Pair[String?, File?]] ome_pairs_input =
    [ ("proteome", prote_ome),
      ("phosphoproteome", phospho_ome),
      ("acetylome", acetyl_ome),
      ("ubiquitylome", ubiquityl_ome) ]
  call select_pairs.panoply_select_all_pairs as ome_pairs { # select extant pairs
    input:
        pairs_input = ome_pairs_input
  }

  Array[Pair[String?, File?]] geneome_pairs_input =
    [ ("rna", "${rna_data}"),
      ("cna", "${cna_data}") ]
  call select_pairs.panoply_select_all_pairs as geneome_pairs { # select extant pairs
    input:
        pairs_input = geneome_pairs_input
  }

  ### NORMALIZE:
  ### Normalize the data first so downstream modules (NMF etc) can run in parallel to main:
  scatter (pair in ome_pairs.pairs) {
    call norm_filt_wdl.panoply_normalize_filter_workflow as norm_filt {
      input:
        input_pome="${pair.right}",
        ome_type="${pair.left}",
        job_identifier="${job_id}-${pair.left}",
        yaml="${yaml}",
        normalizeProteomics=normalizeProteomics,
        filterProteomics=filterProteomics
    }
  }
  
  # Zip Normalization / Filter Output into a labelled Array 
  Array[Pair[String,File]] ome_pairs_norm_filt = zip(norm_filt.output_ome_type , norm_filt.filtered_data_table)


  ### MAIN:
  scatter (pair in ome_pairs_norm_filt) {
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
        input_rna="${rna_data}",
        groups_file="${groups_file}",
        yaml="${yaml}"
    }
  }
  
  # This takes the array of pairs of normalized proteomics data and combines it with the array of pairs of RNA+CNA data for NMF & Blacksheep use:
  Array[Pair[String,File]] all_pairs = flatten([ome_pairs_norm_filt,geneome_pairs.pairs])
  
  ### BLACKSHEEP:
  scatter (pair in all_pairs) {
    call blacksheep_wdl.panoply_blacksheep_workflow as outlier {
      input:
        input_gct = "${pair.right}",
        master_yaml = "${yaml}",
        output_prefix = "${pair.left}",
        type = "${pair.left}",
        groups_file="${if defined(groups_file_blacksheep) then groups_file_blacksheep else groups_file}"
    }
  }

  ### NMF (Multi-omic and Single-omic):
  if ( run_mo_nmf || run_so_nmf ){
    call nmf_wdl.panoply_nmf_workflow as nmf {
      input:
        ome_pairs = all_pairs,

        label = job_id,                     # default parameters & figure colors
        yaml_file = yaml,                   # default parameters & figure colors
        groups_file="${if defined(groups_file_nmf) then groups_file_nmf else groups_file}"

        run_mo_nmf = run_mo_nmf,            # toggle for Multi-omic NMF
        run_so_nmf = run_so_nmf             # toggle for Single-omic NMF
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
          label=job_id,
          groupsFile="${if defined(groups_file_immune) then groups_file_immune else groups_file}"
    }
  }
  
  ## assemble final output combining results from panoply_main, blacksheep immune_analysis and mo_nmf
  call assemble_wdl.panoply_unified_assemble_results {
    input:
      main_full = pome.panoply_full,
      main_summary = pome.summary_and_ssgsea,
      cmap_output = pome.cmap_output,
      cmap_ssgsea_output = pome.cmap_ssgsea_output,
      norm_report = norm_filt.normalize_report,
      rna_corr_report = pome.rna_corr_report,
      cna_corr_report = pome.cna_corr_report,
      omicsev_report = pome.omicsev_report,
      cosmo_report = pome.cosmo_report,
      sampleqc_report = pome.sample_qc_report,
      assoc_report = pome.association_report,
      blacksheep_tar = outlier.blacksheep_tar,
      blacksheep_report = outlier.blacksheep_report,
      nmf_results = nmf.nmf_results,
      nmf_reports = nmf.nmf_reports,
      immune_tar = immune.outputs,
      immune_report = immune.report

  }
  
  output {
    File all_results = panoply_unified_assemble_results.all_results
    File all_reports = panoply_unified_assemble_results.all_reports
  }
 }