#
# Copyright (c) 2025 The Broad Institute, Inc. All rights reserved.
#

import "../panoply_association_workflow/panoply_association_workflow.wdl" as assoc_workflow
import "../panoply_blacksheep_workflow/panoply_blacksheep_workflow.wdl" as blacksheep_wdl
import "../panoply_ssgsea_workflow/panoply_ssgsea_workflow.wdl" as panoply_ssgsea_workflow_wdl
import "../panoply_nmf_internal_workflow/panoply_nmf_internal_workflow.wdl" as nmf_wdl

import "../panoply_immune_analysis_workflow/panoply_immune_analysis_workflow.wdl" as immune_wdl

import "../../tasks/panoply_check_yaml_default/panoply_check_yaml_default.wdl" as check_yaml_default_wdl

workflow panoply_main_internal {

    String job_identifier
    String ome_type
    String? run_ptmsea # "true" or "false"
    String run_nmf = "true"

    ## inputs
    File input_ome
    File yaml

    File groups_file
    File? groups_file_association
    File? groups_file_blacksheep
    File? groups_file_immune
    File? groups_file_nmf

    ## global params
    Float? sample_na_max
    Float? nmiss_factor
    String? duplicate_gene_policy
    String? gene_id_col

    String standalone = "false"
    String geneset_db #this.gseaDB
    String ptm_db #this.ptmseaDB



    #############################
    ### Single-omic Analyses ###
    #############################

    ### Single-Sample GSEA

  call panoply_ssgsea_workflow_wdl.panoply_ssgsea_workflow as ssgsea_ome {
    input:
      preprocess_gct = (ome_type!='rna'), # turn off preprocessing if we have RNA data
      input_ds=input_ome,
      gene_set_database=geneset_db,
      output_prefix=job_identifier,
      level = "gc",
      yaml_file = yaml
  }

  ### Association Analysis
  call assoc_workflow.panoply_association_workflow {
    input: 
      inputData = input_ome, 
      standalone = "true",
      association_groups = select_first([groups_file_association, groups_file]),
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
        label = job_identifier,
        ome_labels=[ome_type],
        ome_gcts=[input_ome],

        yaml_file = yaml,
        groups_file = select_first([groups_file_nmf, groups_file]),
        gene_set_database = geneset_db
    }
  }

  ### BLACKSHEEP:
  call blacksheep_wdl.panoply_blacksheep_workflow as outlier {
    input:
      input_gct = input_ome,
      master_yaml = yaml,
      output_prefix = job_identifier,
      type = ome_type,
      groups_file=select_first([groups_file_blacksheep, groups_file])
  }



  ###################################
  ###   -ome Specific Analyses    ###
  ###################################


  ### IMMUNE: (RNA only)
  if ( ome_type == "rna" ){
    call immune_wdl.panoply_immune_analysis_workflow as immune_analysis {
      input:
          inputData=input_ome,
          standalone="true",
          type=ome_type,
          yaml=yaml,
          analysisDir=job_identifier,
          label=job_identifier,
          groupsFile=select_first([groups_file_immune, groups_file])
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
      call panoply_ssgsea_workflow_wdl.panoply_ssgsea_workflow as ptmsea_ome {
        input:
          preprocess_gct=true,
          input_ds=input_ome,
          gene_set_database=ptm_db,
          output_prefix="${job_identifier}-PTM-SEA",
          level = "ssc",
          yaml_file = yaml
      }
    }
  } 


  output {
    File association_tar = panoply_association_workflow.outputs
    Array[File] ssgsea_assoc_tars = panoply_association_workflow.ssgsea_assoc_tars # association ssgsea results
    File association_report = panoply_association_workflow.report

    File blacksheep_tar = outlier.blacksheep_tar
    File blacksheep_report = outlier.blacksheep_report

    File ssgsea_ome_tar = ssgsea_ome.results
    File ssgsea_ome_report = ssgsea_ome.report

    File? so_nmf_results = so_nmf.nmf_results
    File? so_nmf_figures = so_nmf.nmf_figures
    File? so_nmf_report = so_nmf.nmf_report
    File? so_nmf_ssgsea_tar = so_nmf.nmf_ssgsea_tar
    File? so_nmf_ssgsea_report = so_nmf.nmf_ssgsea_report

    File? immune_analysis_tar = immune_analysis.outputs
    File? immune_analysis_report = immune_analysis.report

    File? ptmsea_tar = ptmsea_ome.results
    File? ptmsea_ome_report = ptmsea_ome.report
  }

}
