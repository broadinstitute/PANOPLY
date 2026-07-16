#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
version 1.0

import "../../tasks/panoply_select_all_pairs/panoply_select_all_pairs.wdl" as select_pairs
import "../panoply_normalize_filter_workflow/panoply_normalize_filter_workflow.wdl" as norm_filt_wdl
import "../panoply_main/panoply_main.wdl" as main_wdl
import "../panoply_main_internal/panoply_main_internal.wdl" as main_internal_wdl
import "../panoply_clumps_ptm_workflow/panoply_clumps_ptm_workflow.wdl" as clumps_wdl
import "../panoply_metaboanalyst_workflow/panoply_metaboanalyst_workflow.wdl" as metab_wdl
import "../panoply_nmf_workflow/panoply_nmf_workflow.wdl" as nmf_wdl
import "../../tasks/panoply_unified_assemble_results/panoply_unified_assemble_results.wdl" as assemble_wdl

import "../../tasks/panoply_check_yaml_default/panoply_check_yaml_default.wdl" as check_yaml_default_wdl


workflow panoply_unified_workflow {
  input {
    File? prote_ome
    File? phospho_ome
    File? acetyl_ome
    File? ubiquityl_ome
    File? nglyco_ome
    File? methyl_ome

    File? metabol_ome

    File? rna_data      #version 1.3 only!
    File? cna_data

    File yaml
    String job_id

    String run_cmap
    Boolean run_mo_nmf #'true' or 'false'
    Boolean run_so_nmf #'true' or 'false'
    String? run_ptmsea
    Boolean? run_clumps
    Boolean? run_metab

    File groups_file
    File? groups_file_nmf
    File? groups_file_metaboanlayst

    File? groups_file_clumpsptm

    File geneset_db
    File ptm_db

    # Normalize specific optional params:
    String? normalizeProteomics # "true" or "false"
    String? filterProteomics # "true" or "false"

    ### Organize omics data into pairs

    # proteomic pairs
  }

  Array[Pair[String?, File?]] ome_pairs_input =
    [ ("proteome", prote_ome),
      ("phosphoproteome", phospho_ome),
      ("acetylome", acetyl_ome),
      ("ubiquitylome", ubiquityl_ome),
      ("nglycoproteome", nglyco_ome),
      ("methylation", methyl_ome) ]

  call select_pairs.panoply_select_all_pairs as ome_pairs { # select extant pairs
    input:
        pairs_input = ome_pairs_input
  }

  # genomic pairs
  Array[Pair[String?, File?]] genome_pairs_input =
    [ ("rna", rna_data),
      ("cna", cna_data) ]
  call select_pairs.panoply_select_all_pairs as genome_pairs { # select extant pairs
    input:
        pairs_input = genome_pairs_input
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
  
  # This takes the array of pairs of normalized proteomics data and combines it with the array of pairs of RNA+CNA data for NMF & Blacksheep use:
  Array[Pair[String,File]] all_pairs = flatten([ome_pairs_norm_filt,genome_pairs.pairs])


  ### MAIN:
  scatter (pair in ome_pairs_norm_filt) {
    call main_wdl.panoply_main as pome {
      input:
        ## include all required arguments from above
        input_pome=pair.right,
        ome_type=pair.left,
        job_identifier="${job_id}-${pair.left}",
        geneset_db=geneset_db,
        run_ptmsea="${run_ptmsea}",
        ptm_db=ptm_db,
        run_cmap = "${run_cmap}",
        run_nmf = "false",
        input_cna=cna_data,
        input_rna=rna_data,
        groups_file=groups_file,
        yaml=yaml
    }
  }

  ### MAIN for RNA:
  if (defined(rna_data)) {
    call main_internal_wdl.panoply_main_internal as rna {
      input:
        input_ome=select_first([rna_data]),
        ome_type="rna",
        job_identifier="${job_id}-rna",
        geneset_db=geneset_db,
        run_ptmsea=false,
        ptm_db=ptm_db,
        run_nmf = "false",
        groups_file=groups_file,
        yaml=yaml
    }
  }

  ### ClumpsPTM
  # check yaml default for run.clumpsptm (Terra param takes precedence)
  call check_yaml_default_wdl.panoply_check_yaml_default as check_clumpsptm_default {
    input:
      param = run_clumps,
      yaml = yaml,
      param_lookup = "run.clumpsptm"
  }
  if ( check_clumpsptm_default.param_boolean ){
    call clumps_wdl.panoply_clumps_ptm_workflow as clumps_ptm {
      input:
        pSTY_gct = phospho_ome,
        acK_gct = acetyl_ome,
        ubK_gct = ubiquityl_ome,
        groupsFile = select_first([groups_file_clumpsptm]), # NOTE: this will fail somewhat lazily if groups_file_clumpsptm is not provided
        output_prefix = job_id,
        yaml_file = yaml
    }
  }

  ### MetaboAnalyst
  # check yaml default for run.metab (Terra param takes precedence)
  call check_yaml_default_wdl.panoply_check_yaml_default as check_metab_default {
    input:
      param = run_metab,
      yaml = yaml,
      param_lookup = "run.metab"
  }
  if ( defined(metabol_ome) && check_metab_default.param_boolean ){
    # Proteome and Transcriptome pair
    Array[Pair[String?, File?]] pg_pairs_input =
      [ ("proteome", prote_ome),
        ("rna", rna_data) ]
    call select_pairs.panoply_select_all_pairs as pg_pairs { # select extant pairs
      input:
          pairs_input = pg_pairs_input
    }

    scatter (pair in pg_pairs.pairs) {
      call metab_wdl.panoply_metaboanalyst_workflow as metab {
        input:
          meta_gct = select_first([metabol_ome]),
          omic_gct = pair.right,
          ome_type = pair.left,
          output_prefix = "${job_id}-${pair.left}",
          groups_file = select_first([groups_file_metaboanlayst, groups_file]),
          yaml_file = yaml
      }
    }
  }

  ### NMF (Multi-omic and Single-omic):
  if ( run_mo_nmf || run_so_nmf ){
    call nmf_wdl.panoply_nmf_workflow as nmf {
      input:
        ome_pairs = all_pairs,

        label = job_id,                     # default parameters & figure colors
        yaml_file = yaml,                   # default parameters & figure colors
        groups_file=select_first([groups_file_nmf, groups_file]),

        gene_set_database=geneset_db,

        run_mo_nmf = run_mo_nmf,            # toggle for Multi-omic NMF
        run_so_nmf = run_so_nmf             # toggle for Single-omic NMF
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
      ssgsea_rna_report = rna.ssgsea_ome_report,
      ssgsea_ome_report = pome.ssgsea_ome_report,
      omicsev_report = pome.omicsev_report,
      cosmo_report = pome.cosmo_report,
      sampleqc_report = pome.sample_qc_report,
      assoc_report = pome.association_report,
      ptmsea_ome_report = pome.ptmsea_ome_report,
      blacksheep_report = pome.blacksheep_report,
      clumpsptm_results = clumps_ptm.clumps_ptm_results,
      clumpsptm_report = clumps_ptm.clumps_ptm_report,
      metaboanalyst_results = metab.metaboanalyst_tar,
      metaboanalyst_reports = metab.metaboanalyst_report,
      nmf_results = nmf.nmf_results,
      nmf_reports = nmf.nmf_reports,
      immune_report = rna.immune_analysis_report,
      rna_blacksheep_report = rna.blacksheep_report,
      rna_assoc_report = rna.association_report
  }
  
  output {
    File all_results = panoply_unified_assemble_results.all_results
    File all_reports = panoply_unified_assemble_results.all_reports
  }
 }