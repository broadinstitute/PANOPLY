#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#

import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_rna_protein_correlation/versions/4/plain-WDL/descriptor" as rna_prot_corr_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_harmonize/versions/4/plain-WDL/descriptor" as harmonize_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_sampleqc/versions/4/plain-WDL/descriptor" as sampleqc_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_cna_setup/versions/4/plain-WDL/descriptor" as cna_setup_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_cna_correlation/versions/2/plain-WDL/descriptor" as cna_corr_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_association/versions/5/plain-WDL/descriptor" as assoc_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_accumulate/versions/4/plain-WDL/descriptor" as accum_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_download/versions/3/plain-WDL/descriptor" as download_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_ssgsea/versions/14/plain-WDL/descriptor" as ssgsea_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_rna_protein_correlation_report/versions/3/plain-WDL/descriptor" as rna_corr_report_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_cna_correlation_report/versions/4/plain-WDL/descriptor" as cna_corr_report_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_sampleqc_report/versions/4/plain-WDL/descriptor" as sampleqc_report_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_association_report/versions/6/plain-WDL/descriptor" as assoc_report_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_cmap_analysis/versions/5/plain-WDL/descriptor" as cmap_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_omicsev/versions/1/plain-WDL/descriptor" as omicsev_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_so_nmf_gct/versions/9/plain-WDL/descriptor" as so_nmf_wdl


workflow panoply_main {

  String job_identifier
  String ome_type
  String run_ptmsea # "true" or "false"
  File sample_annotation
  String run_cmap   # "true" or "false"
  String? run_nmf = "true"

  ## inputs
  File input_pome
  File input_rna_v3
  File input_cna
  File yaml

  File? cna_groups
  File? association_groups

  ## cmap inputs
  Int cmap_n_permutations = 10
  File? cmap_enrichment_groups
  File subset_list_file = "gs://fc-de501ca1-0ae7-4270-ae76-6c99ea9a6d5b/cmap-data/cmap-data-subsets-index.txt"
  File cmap_level5_data = "gs://fc-de501ca1-0ae7-4270-ae76-6c99ea9a6d5b/cmap-data/annotated_GSE92742_Broad_LINCS_Level5_COMPZ_geneKDsubset_n36720x12328.gctx"
  File? annotation_pathway_db #this.gseaDB
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
  
  call rna_prot_corr_wdl.panoply_rna_protein_correlation {
    input:
      inputData = input_pome,
      type = ome_type,
      rnaExpr = input_rna_v3,
      analysisDir = job_identifier,
      standalone = "true",
      yaml = yaml
  }

  call ssgsea_wdl.panoply_ssgsea as ssgsea_rna {
    input:
      input_ds = input_rna_v3,
      gene_set_database = geneset_db,
      output_prefix = job_identifier,
      level = "gc",
      yaml_file = yaml
  }

  call ssgsea_wdl.panoply_ssgsea as ssgsea_ome {
    input:
      input_ds = input_pome,
      gene_set_database = geneset_db,
      output_prefix = job_identifier,
      level = "gc",
      yaml_file = yaml
  }
  
  if ( run_ptmsea == "true" ){
    if ( ome_type == "phosphoproteome" ){
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

  call rna_corr_report_wdl.panoply_rna_protein_correlation_report {
    input:
      tarball = panoply_rna_protein_correlation.outputs,
      config_yaml = yaml,
      label = job_identifier,
      type = ome_type,
      tmpDir = "tmp"
  }

  call harmonize_wdl.panoply_harmonize {
    input:
      inputData = panoply_rna_protein_correlation.outputs,
      rnaExpr = input_rna_v3,
      cnaExpr = input_cna,
      standalone = standalone,
      type = ome_type,
      yaml = yaml,
      na_max=na_max,
      duplicate_gene_policy=duplicate_gene_policy,
      gene_id_col=gene_id_col
  }
  
  call omicsev_wdl.panoply_omicsev {
    input:
      yaml_file = yaml,
      STANDALONE = standalone,
      do_function_prediction = false,
      panoply_harmonize_tar_file = panoply_harmonize.outputs,
      label = job_identifier
  }

  call sampleqc_wdl.panoply_sampleqc {
    input:
      tarball = panoply_omicsev.outputs,
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

  call cna_setup_wdl.panoply_cna_setup {
    input:
      tarball = panoply_sampleqc.outputs,
      groupsFile = cna_groups,
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
  
  call assoc_wdl.panoply_association {
    input: 
      inputData = panoply_cna_correlation.outputs, 
      groupsFile = association_groups,
      type = ome_type,
      standalone = standalone,
      yaml = yaml,
      sample_na_max=sample_na_max,
      nmiss_factor=nmiss_factor,
      duplicate_gene_policy=duplicate_gene_policy,
      gene_id_col=gene_id_col
  }

  call accum_wdl.panoply_accumulate as accumulate_assoc {
    input:
      input_tar = panoply_association.outputs,
      module = "association"
  } 

  Array[File] list_gct_assoc = accumulate_assoc.list_gct
  scatter (f in list_gct_assoc){
    call ssgsea_wdl.panoply_ssgsea as ssgsea_assoc {
      input:
        input_ds = "${f}",
        gene_set_database = geneset_db,
        output_prefix = job_identifier,
        level = "gc",
        yaml_file = yaml
    }
  }
  
  if ( run_cmap == "true" ){
    if ( ome_type == "proteome" ) {
      call cmap_wdl.run_cmap_analysis {
        input:
          CNAcorr_tarball = panoply_cna_correlation.outputs,
          subset_list_file = subset_list_file,
          cmap_level5_data = cmap_level5_data,
          annotation_pathway_db = annotation_pathway_db, 
          subset_bucket = subset_bucket,
          n_permutations = cmap_n_permutations,
            cmap_enrichment_groups = cmap_enrichment_groups,
            yaml = yaml
        
      }
    }
  }

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

  call download_wdl.panoply_download {
    input:
      association_tar = panoply_association.outputs,
      ssgsea_ome_tar = ssgsea_ome.results,
      ssgsea_rna_tar = ssgsea_rna.results,
      analysisDir = job_identifier,
      ssgsea_assoc_tars = ssgsea_assoc.results,
      ptmsea = ptmsea_ome.results,
      so_nmf_tar = so_nmf.nmf_clust,
      so_nmf_ssgsea_tar = so_nmf.nmf_ssgsea,
      output_prefix = ome_type
  }
  
  call assoc_report_wdl.panoply_association_report {
    input:
      input_tar = panoply_download.full,
      master_yaml = yaml,
      label = "${job_identifier}-full-results",
      type = ome_type
  }

  output {
    File summary_and_ssgsea = panoply_download.summary
    File panoply_full = panoply_download.full
    File rna_corr_report = panoply_rna_protein_correlation_report.report
    File cna_corr_report = panoply_cna_correlation_report.report
    File omicsev_report = panoply_omicsev.report
    File sample_qc_report = panoply_sampleqc_report.report
    File association_report = panoply_association_report.report_out
    File? so_nmf_report = so_nmf.nmf_clust_report
    File? so_nmf_ssgsea_report = so_nmf.nmf_ssgsea_report
    File? cmap_output = run_cmap_analysis.outputs
    File? cmap_ssgsea_output = run_cmap_analysis.ssgseaOutput
  }

}