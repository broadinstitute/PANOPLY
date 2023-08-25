#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#

# association analysis
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_association_workflow/versions/1/plain-WDL/descriptor" as assoc_workflow_wdl

# ssgsea and ptmsea
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_ssgsea/versions/6/plain-WDL/descriptor" as ssgsea_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_ssgsea_report/versions/8/plain-WDL/descriptor" as ssgsea_report_wdl

# nmf
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_so_nmf_gct/versions/17/plain-WDL/descriptor" as so_nmf_wdl

# misc
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_check_yaml_default/versions/2/plain-WDL/descriptor" as check_yaml_default_wdl
#import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_omicsev/versions/21/plain-WDL/descriptor" as omicsev_wdl
#import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_cosmo/versions/11/plain-WDL/descriptor" as cosmo_wdl


workflow panoply_main_NG {

  String job_identifier
  String ome_type
  String? run_ptmsea # "true" or "false"
  File sample_annotation
  String? run_nmf = "true"

  ## inputs
  File input_pome
  File yaml

  File association_groups
  
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

  String standalone = "true"
  String geneset_db #this.gseaDB
  String ptm_db #this.ptmseaDB
  


  ### Association Analysis
  call assoc_workflow_wdl.panoply_association_workflow as assoc {
    input: 
      input_pome = input_pome, 
      association_groups = association_groups,
      geneset_db=geneset_db,
      ome_type = ome_type,
      job_identifier = job_identifier,
      yaml = yaml,
      sample_na_max=sample_na_max,
      nmiss_factor=nmiss_factor,
      duplicate_gene_policy=duplicate_gene_policy,
      gene_id_col=gene_id_col
  }



  ### SSGSEA

  call ssgsea_wdl.panoply_ssgsea as ssgsea_ome {
    input:
      input_ds = input_pome,
      gene_set_database = geneset_db,
      output_prefix = job_identifier,
      level = "gc",
      yaml_file = yaml
  }
  call ssgsea_report_wdl.panoply_ssgsea_report as ssgsea_ome_report {
    input:
      tarball = ssgsea_ome.results,
      cfg_yaml = yaml,
      label = job_identifier
  }


  ### PTMSEA

  if ( ome_type == "phosphoproteome" ){

    # check yaml default for run.ptmsea (Terra param takes precedence)
    call check_yaml_default_wdl.panoply_check_yaml_default as check_ptmsea_default {
      input:
        param = run_ptmsea,
        yaml = yaml,
        param_lookup = "run.ptmsea"
    }

    if ( check_ptmsea_default.param_boolean ){
      call ssgsea_wdl.panoply_ssgsea as ptmsea {
        input:
          input_ds = input_pome,
          gene_set_database = ptm_db,
          output_prefix = job_identifier,
          level = "ssc",
          yaml_file = yaml
      }
      call ssgsea_report_wdl.panoply_ssgsea_report as ptmsea_report {
        input:
          tarball = ptmsea.results,
          cfg_yaml = yaml,
          label = job_identifier,
          is_ptmsigdb = true # set ptmsea to true to use ptmsigdb signatures in heatmaps
      }
    }
  }



  ### SO-NMF Clustering 

  if ( run_nmf == "true" ){
    call so_nmf_wdl.panoply_so_nmf_gct_workflow as so_nmf {
      input:
      yaml_file = yaml,
      label = job_identifier,
      ome = input_pome,
      ome_type = ome_type,
      gene_set_database = geneset_db
    }
  }




#  call cosmo_wdl.panoply_cosmo_workflow {
#    input:
#      STANDALONE = standalone,
#      yaml_file = yaml,
#      panoply_harmonize_tar = panoply_harmonize.outputs,
#      label = job_identifier,
#      ome_type = ome_type
#  }
  
#  call omicsev_wdl.panoply_omicsev {
#    input:
#      yaml_file = yaml,
#      STANDALONE = standalone,
#      do_function_prediction = false,
#      panoply_harmonize_tar_file = panoply_harmonize.outputs,
#      label = job_identifier,
#      ome_type = ome_type
#  }

  output {
    File association_contrasts = assoc.contrasts
    File association_report = assoc.report
    
    File ssgsea_results = ssgsea_ome.results
    File ssgsea_report = ssgsea_ome_report.report
    
    File? ptmsea_results = ptmsea.results
    File? ptmsea_report = ptmsea_report.report

    File? so_nmf_results = so_nmf.nmf_clust
    File? so_nmf_report = so_nmf.nmf_clust_report
    File? so_nmf_ssgsea_report = so_nmf.nmf_ssgsea_report
    
    #File omicsev_report = panoply_omicsev.report
    #File? cosmo_report = panoply_cosmo_workflow.cosmo_report
  }

}