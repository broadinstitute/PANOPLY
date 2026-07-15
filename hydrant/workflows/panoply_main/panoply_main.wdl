#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
version 1.0

import "../panoply_main_internal/panoply_main_internal.wdl" as panoply_main_internal
## Proteogenomic
import "../../tasks/panoply_harmonize/panoply_harmonize.wdl" as harmonize_wdl
import "../../tasks/panoply_rna_protein_correlation/panoply_rna_protein_correlation.wdl" as rna_prot_corr_wdl
import "../../tasks/panoply_rna_protein_correlation_report/panoply_rna_protein_correlation_report.wdl" as rna_corr_report_wdl
import "../../tasks/panoply_cna_setup/panoply_cna_setup.wdl" as cna_setup_wdl
import "../../tasks/panoply_cna_correlation/panoply_cna_correlation.wdl" as cna_corr_wdl
import "../../tasks/panoply_cna_correlation_report/panoply_cna_correlation_report.wdl" as cna_corr_report_wdl
import "../../tasks/panoply_cmap_analysis/panoply_cmap_analysis.wdl" as cmap_wdl
## Sample-QC
import "../../tasks/panoply_sampleqc/panoply_sampleqc.wdl" as sampleqc_wdl
import "../../tasks/panoply_sampleqc_report/panoply_sampleqc_report.wdl" as sampleqc_report_wdl
import "../../tasks/panoply_cosmo/panoply_cosmo.wdl" as cosmo_wdl
import "../../tasks/panoply_omicsev/panoply_omicsev.wdl" as omicsev_wdl
## Support
import "../../tasks/panoply_download/panoply_download.wdl" as download_wdl


workflow panoply_main {

  input {
    String job_identifier
    String ome_type
    String? run_ptmsea # "true" or "false"
    String run_cmap   # "true" or "false"
    String run_nmf = "true"
    String run_omicsev = "true"

    ## inputs
    File input_pome
    File? input_rna
    File? input_cna
    File yaml

    File groups_file
    File? groups_file_association
    File? groups_file_blacksheep
    File? groups_file_cmap_enrichment
    File? groups_file_immune
    File? groups_file_nmf

    File? cna_corr_groupsFile # DO NOT use groups_file by default

    ## cmap inputs
    Int cmap_n_permutations = 10
    File subset_list_file = "gs://fc-de501ca1-0ae7-4270-ae76-6c99ea9a6d5b/cmap-data/cmap-data-subsets-index.txt"
    File cmap_level5_data = "gs://fc-de501ca1-0ae7-4270-ae76-6c99ea9a6d5b/cmap-data/annotated_GSE92742_Broad_LINCS_Level5_COMPZ_geneKDsubset_n36720x12328.gctx"
    File? geneset_db_cmap # optional override for CMAP geneset_db
    String subset_bucket = "gs://fc-de501ca1-0ae7-4270-ae76-6c99ea9a6d5b/cmap-data/cmap-data-subsets"

    ## global params
    Float? na_max
    Float? sample_na_max
    Float? nmiss_factor
    String? duplicate_gene_policy
    String? gene_id_col

    String standalone = "false"
    String geneset_db #this.gseaDB
    String ptm_db #this.ptmseaDB


    ###################################
    ###   Proteogenomics Analyses   ###
    ###################################
  }

  if (defined(input_rna) && defined(input_cna)) {

    ### RNA Correlation
    call rna_prot_corr_wdl.panoply_rna_protein_correlation {
      input:
        inputData = input_pome,
        type = ome_type,
        rnaExpr = input_rna,
        analysisDir = job_identifier,
        standalone = "true",
        yaml = yaml
    }

    call rna_corr_report_wdl.panoply_rna_protein_correlation_report {
      input:
        tarball = panoply_rna_protein_correlation.outputs,
        config_yaml = yaml,
        label = job_identifier,
        type = ome_type,
        tmpDir = "tmp"
    }


    ### Harmonize Datasets
    call harmonize_wdl.panoply_harmonize {
      input:
        inputData = panoply_rna_protein_correlation.outputs,
        rnaExpr = input_rna,
        cnaExpr = input_cna,
        standalone = standalone,
        type = ome_type,
        yaml = yaml,
        na_max=na_max,
        duplicate_gene_policy=duplicate_gene_policy,
        gene_id_col=gene_id_col
    }
    
    ### COSMO
    call cosmo_wdl.panoply_cosmo_workflow {
      input:
        STANDALONE = standalone,
        yaml_file = yaml,
        panoply_harmonize_tar = panoply_harmonize.outputs,
        label = job_identifier,
        ome_type = ome_type
    }
    
    ### Omics EV
    if ( run_omicsev == "true" ){ 
      call omicsev_wdl.panoply_omicsev {
        input:
          yaml_file = yaml,
          STANDALONE = standalone,
          do_function_prediction = false,
          panoply_harmonize_tar_file = panoply_harmonize.outputs,
          label = job_identifier,
          ome_type = ome_type
      }
    }


    ### Sample QC
    call sampleqc_wdl.panoply_sampleqc {
      input:
        tarball = panoply_harmonize.outputs,
        type = ome_type,
        yaml = yaml
    }

    call sampleqc_report_wdl.panoply_sampleqc_report {
      input:
        tarball = panoply_sampleqc.outputs,
        type = ome_type,
        label = job_identifier,
        tmpDir = "tmp"
    }

    ### CNA Correlation
    call cna_setup_wdl.panoply_cna_setup {
      input:
        tarball = panoply_sampleqc.outputs,
        groupsFile = cna_corr_groupsFile,
        type = ome_type,
        yaml = yaml
    }

    call cna_corr_wdl.panoply_cna_correlation {
      input:
        tarball = panoply_cna_setup.outputs,
        type = ome_type,
        yaml = yaml
    }

    call cna_corr_report_wdl.panoply_cna_correlation_report {
      input:
        tarball = panoply_cna_correlation.outputs,
        config_yaml = yaml,
        type = ome_type,
        label = job_identifier,
        tmpDir = "tmp"
    }

    ### CMAP Analysis (proteome only)
    if ( run_cmap == "true" ){
      if ( ome_type == "proteome" ) {
        call cmap_wdl.run_cmap_analysis {
          input:
            CNAcorr_tarball = panoply_cna_correlation.outputs,
            subset_list_file = subset_list_file,
            cmap_level5_data = cmap_level5_data,
            annotation_pathway_db = select_first([geneset_db_cmap, geneset_db]),
            subset_bucket = subset_bucket,
            n_permutations = cmap_n_permutations,
            cmap_enrichment_groups = select_first([groups_file_cmap_enrichment, groups_file]),
            yaml = yaml
          
        }
      }
    }
  }

  #############################
  ### Non-Genomics Analyses ###
  #############################

  call panoply_main_internal.panoply_main_internal as main_internal {
    input:
      job_identifier = job_identifier,
      ome_type = ome_type,
      run_ptmsea = run_ptmsea,
      run_nmf = run_nmf,
      input_ome = input_pome,
      yaml = yaml,
      groups_file = groups_file,
      groups_file_association = groups_file_association,
      groups_file_blacksheep = groups_file_blacksheep,
      groups_file_immune = groups_file_immune,
      groups_file_nmf = groups_file_nmf,
      sample_na_max = sample_na_max,
      nmiss_factor = nmiss_factor,
      duplicate_gene_policy = duplicate_gene_policy,
      gene_id_col = gene_id_col,
      standalone = standalone,
      geneset_db = geneset_db,
      ptm_db = ptm_db
  }



  #############################
  ###   Compiile Results    ###
  #############################
  call download_wdl.panoply_download {
    input:
      association_tar = main_internal.association_tar,
      ssgsea_assoc_tars = main_internal.ssgsea_assoc_tars, # association ssgsea results
      ssgsea_ome_tar = main_internal.ssgsea_ome_tar,
      blacksheep_tar = main_internal.blacksheep_tar,
      so_nmf_results = main_internal.so_nmf_results,
      so_nmf_figures = main_internal.so_nmf_figures,
      so_nmf_ssgsea_tar = main_internal.so_nmf_ssgsea_tar,
      immune_analysis_tar = main_internal.immune_analysis_tar,
      ptmsea = main_internal.ptmsea_tar,
      omicsev_tar = panoply_omicsev.outputs,
      cosmo_tar = panoply_cosmo_workflow.cosmo_tar,
      cna_corr_tar = panoply_cna_correlation.outputs, # contains all results in non-standalone
      analysisDir = job_identifier,
      output_prefix = ome_type
  }

  output {
    File summary_and_ssgsea = panoply_download.summary
    File panoply_full = panoply_download.full

    # multiomic analyses
    File? cmap_output = run_cmap_analysis.outputs
    File? cmap_ssgsea_output = run_cmap_analysis.ssgseaOutput
    File? rna_corr_report = panoply_rna_protein_correlation_report.report
    File? cna_corr_report = panoply_cna_correlation_report.report
    File? omicsev_report = panoply_omicsev.report
    File? cosmo_report = panoply_cosmo_workflow.cosmo_report
    File? sample_qc_report = panoply_sampleqc_report.report

    # single-omic analyses
    File association_report = main_internal.association_report
    File blacksheep_report = main_internal.blacksheep_report
    File ssgsea_ome_report = main_internal.ssgsea_ome_report
    File? so_nmf_report = main_internal.so_nmf_report
    File? so_nmf_ssgsea_report = main_internal.so_nmf_ssgsea_report
    # -ome specific analyses
    File? immune_report = main_internal.immune_analysis_report
    File? ptmsea_ome_report = main_internal.ptmsea_ome_report
  }

}