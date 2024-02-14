#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#

import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_rna_protein_correlation/versions/18/plain-WDL/descriptor" as rna_prot_corr_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_harmonize/versions/17/plain-WDL/descriptor" as harmonize_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_sampleqc/versions/17/plain-WDL/descriptor" as sampleqc_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_cna_setup/versions/18/plain-WDL/descriptor" as cna_setup_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_cna_correlation/versions/13/plain-WDL/descriptor" as cna_corr_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_rna_protein_correlation_report/versions/16/plain-WDL/descriptor" as rna_corr_report_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_cna_correlation_report/versions/18/plain-WDL/descriptor" as cna_corr_report_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_sampleqc_report/versions/17/plain-WDL/descriptor" as sampleqc_report_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_cmap_analysis/versions/17/plain-WDL/descriptor" as cmap_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_check_yaml_default/versions/28/plain-WDL/descriptor" as check_yaml_default_wdl

import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_association_workflow/versions/19/plain-WDL/descriptor" as assoc_workflow
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_ssgsea/versions/37/plain-WDL/descriptor" as ssgsea_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_omicsev/versions/4/plain-WDL/descriptor" as omicsev_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_nmf_internal_workflow/versions/16/plain-WDL/descriptor" as nmf_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_download/versions/17/plain-WDL/descriptor" as download_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_cosmo/versions/3/plain-WDL/descriptor" as cosmo_wdl


workflow panoply_main {

  String job_identifier
  String ome_type
  String? run_ptmsea # "true" or "false"
  String run_cmap   # "true" or "false"
  String? run_nmf = "true"

  ## inputs
  File input_pome
  File? input_rna
  File? input_cna
  File yaml

  File groups_file
  File? groups_file_association
  File? groups_file_cmap_enrichment
  File? groups_file_nmf

  File? cna_corr_groupsFile # DO NOT use groups_file by default

  ## cmap inputs
  Int cmap_n_permutations = 10
  File subset_list_file = "gs://fc-de501ca1-0ae7-4270-ae76-6c99ea9a6d5b/cmap-data/cmap-data-subsets-index.txt"
  File cmap_level5_data = "gs://fc-de501ca1-0ae7-4270-ae76-6c99ea9a6d5b/cmap-data/annotated_GSE92742_Broad_LINCS_Level5_COMPZ_geneKDsubset_n36720x12328.gctx"
  File? geneset_db_cmap # optional override for CMAP geneset_db
  String subset_bucket = "gs://fc-de501ca1-0ae7-4270-ae76-6c99ea9a6d5b/cmap-data/cmap-data-subsets"
  
  ## global params
  Int? ndigits
  Float? na_max
  Float? sample_na_max
  Float? min_numratio_fraction
  Float? nmiss_factor
  Float? sd_filter_threshold
  String? duplicate_gene_policy
  String? gene_id_col
  String? organism

  String standalone = "false"
  String geneset_db #this.gseaDB
  String ptm_db #this.ptmseaDB
  

  #############################
  ###   Genomics Analyses   ###
  #############################

  if (defined(input_rna) && defined(input_cna)) {

    ### Single-Sample GSEA (on RNA)
    call ssgsea_wdl.panoply_ssgsea as ssgsea_rna {
      input:
        input_ds = input_rna,
        gene_set_database = geneset_db,
        output_prefix = job_identifier,
        level = "gc",
        yaml_file = yaml
    }

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
    call omicsev_wdl.panoply_omicsev {
      input:
        yaml_file = yaml,
        STANDALONE = standalone,
        do_function_prediction = false,
        panoply_harmonize_tar_file = panoply_harmonize.outputs,
        label = job_identifier,
        ome_type = ome_type
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
            annotation_pathway_db = "${if defined(geneset_db_cmap) then geneset_db_cmap else geneset_db}",
            subset_bucket = subset_bucket,
            n_permutations = cmap_n_permutations,
            cmap_enrichment_groups = "${if defined(groups_file_cmap_enrichment) then groups_file_cmap_enrichment else groups_file}",
            yaml = yaml
          
        }
      }
    }
  }

  #############################
  ### Non-Genomics Analyses ###
  #############################

  ### Single-Sample GSEA (on pome)
  call ssgsea_wdl.panoply_ssgsea as ssgsea_ome {
    input:
      input_ds = input_pome,
      gene_set_database = geneset_db,
      output_prefix = job_identifier,
      level = "gc",
      yaml_file = yaml
  }

  ### Association Analysis
  call assoc_workflow.panoply_association_workflow {
    input: 
      inputData = input_pome, 
      standalone = "true",
      association_groups = "${if defined(groups_file_association) then groups_file_association else groups_file}",
      geneset_db=geneset_db,
      ome_type = ome_type,
      job_identifier = job_identifier,
      yaml = yaml,
      sample_na_max=sample_na_max,
      nmiss_factor=nmiss_factor,
      duplicate_gene_policy=duplicate_gene_policy,
      gene_id_col=gene_id_col
  }

  ### NMF Analysis
  if ( run_nmf == "true" ){
    call nmf_wdl.panoply_nmf_internal_workflow as so_nmf {
      input:
        label = "${job_identifier}_${ome_type}",
        ome_labels=[ome_type],
        ome_gcts=[input_pome],

        yaml_file = yaml,
        groups_file = "${if defined(groups_file_nmf) then groups_file_nmf else groups_file}",
        gene_set_database = geneset_db
    }
  }

  ### PTMSEA (phosphoproteome only)
  if ( ome_type == "phosphoproteome" ){

    # check yaml default for run.ptmsea (Terra param takes precedence)
    call check_yaml_default_wdl.panoply_check_yaml_default as check_ptmsea_default {
      input:
        param = run_ptmsea,
        yaml = yaml,
        param_lookup = "run.ptmsea"
    }

    if ( check_ptmsea_default.param_boolean ){
      call ssgsea_wdl.panoply_ssgsea as ptmsea_ome {
        input:
          input_ds = input_pome,
          gene_set_database = ptm_db,
          output_prefix = job_identifier,
          level = "ssc",
          yaml_file = yaml
      }
    }
  } 



  #############################
  ###   Compiile Results    ###
  #############################
  call download_wdl.panoply_download {
    input:
      association_tar = panoply_association_workflow.outputs,
      cna_corr_tar = panoply_cna_correlation.outputs, # contains all results in non-standalone
      ssgsea_ome_tar = ssgsea_ome.results,
      ssgsea_rna_tar = ssgsea_rna.results,
      omicsev_tar = panoply_omicsev.outputs,
      cosmo_tar = panoply_cosmo_workflow.cosmo_tar,
      analysisDir = job_identifier,
      ssgsea_assoc_tars = panoply_association_workflow.ssgsea_assoc_tars, # association ssgsea results
      ptmsea = ptmsea_ome.results,
      so_nmf_results = so_nmf.nmf_results,
      so_nmf_figures = so_nmf.nmf_figures,
      so_nmf_ssgsea_tar = so_nmf.nmf_ssgsea_tar,
      output_prefix = ome_type
  }

  output {
    File summary_and_ssgsea = panoply_download.summary
    File panoply_full = panoply_download.full
    File? rna_corr_report = panoply_rna_protein_correlation_report.report
    File? cna_corr_report = panoply_cna_correlation_report.report
    File? omicsev_report = panoply_omicsev.report
    File? cosmo_report = panoply_cosmo_workflow.cosmo_report
    File? sample_qc_report = panoply_sampleqc_report.report
    File association_report = panoply_association_workflow.report
    File? so_nmf_report = so_nmf.nmf_report
    File? so_nmf_ssgsea_report = so_nmf.nmf_ssgsea_report
    File? cmap_output = run_cmap_analysis.outputs
    File? cmap_ssgsea_output = run_cmap_analysis.ssgseaOutput
  }

}