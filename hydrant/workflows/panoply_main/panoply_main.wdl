import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_rna_protein_correlation/versions/1/plain-WDL/descriptor" as rna_prot_corr_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_harmonize/versions/1/plain-WDL/descriptor" as harmonize_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_sampleqc/versions/1/plain-WDL/descriptor" as sampleqc_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_cna_setup/versions/1/plain-WDL/descriptor" as cna_setup_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_cna_correlation/versions/1/plain-WDL/descriptor" as cna_corr_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_association/versions/1/plain-WDL/descriptor" as assoc_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_accumulate/versions/1/plain-WDL/descriptor" as accum_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_download/versions/1/plain-WDL/descriptor" as download_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:pgdac_ssgsea/versions/8/plain-WDL/descriptor" as ssgsea_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:pgdac_rna_protein_correlation_report/versions/1/plain-WDL/descriptor" as rna_corr_report_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:pgdac_cna_correlation_report/versions/1/plain-WDL/descriptor" as cna_corr_report_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:pgdac_sampleqc_report/versions/1/plain-WDL/descriptor" as sampleqc_report_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:pgdac_cons_clust/versions/3/plain-WDL/descriptor" as cons_clust_wdl


workflow panoply_main {

  String job_identifier
  String ome_type
  String run_ptmsea
  File sample_annotation

  ## inputs
  File input_pome
  File input_rna
  File input_rna_v3
  File input_cna

  ## FDRs
  Float fdr_for_reports
  Float fdr_for_cna_corr

  String? data_sub_type
  File? add_parameters
  File? cna_groups
  File? association_groups
  File? cluster_enrichment_groups

  String standalone = "false"
  String geneset_db = "gs://fc-e9c1f751-c433-464e-936e-c795faf4eca0/msigdb_v7.0_h.all.v7.0.symbols.gmt"
  String ptm_db = "gs://fc-e9c1f751-c433-464e-936e-c795faf4eca0/ptm.sig.db.all.uniprot.human.v1.9.0.gmt"

  call rna_prot_corr_wdl.panoply_rna_protein_correlation {
    input:
      inputData = input_pome,
      type = ome_type,
      rnaExpr = input_rna,
      analysisDir = job_identifier,
      standalone = "true",
      params = add_parameters
  }

  call ssgsea_wdl.pgdac_ssgsea as ssgsea_rna {
    input:
      input_ds = input_rna_v3,
      gene_set_database = geneset_db,
      output_prefix = job_identifier,
      level = "gc"
  }

  call ssgsea_wdl.pgdac_ssgsea as ssgsea_ome {
    input:
      input_ds = input_pome,
      gene_set_database = geneset_db,
      output_prefix = job_identifier,
      level = "gc"
  }
  
  if ( run_ptmsea == "true" ){
    if ( ome_type == "phosphoproteome" ){
      call ssgsea_wdl.pgdac_ssgsea as ptmsea_ome {
        input:
          input_ds = input_pome,
          gene_set_database = ptm_db,
          output_prefix = job_identifier,
          level = "ssc"
      }
    }
  } 

  call rna_corr_report_wdl.pgdac_rna_protein_correlation_report {
    input:
      tarball = panoply_rna_protein_correlation.outputs,
      label = job_identifier,
      fdr = fdr_for_reports,
      type = ome_type,
      tmpDir = "tmp"
  }

  call harmonize_wdl.panoply_harmonize {
    input:
      inputData = panoply_rna_protein_correlation.outputs,
      rnaExpr = input_rna,
      cnaExpr = input_cna,
      standalone = standalone,
      type = ome_type,
      subType = data_sub_type,
      params = add_parameters
  }

  call sampleqc_wdl.panoply_sampleqc {
    input:
      tarball = panoply_harmonize.outputs,
      type = ome_type,
      subType = data_sub_type,
      params = add_parameters
  }

  call sampleqc_report_wdl.pgdac_sampleqc_report {
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
      subType = data_sub_type,
      params = add_parameters
  }

  call cna_corr_wdl.panoply_cna_correlation {
    input:
      tarball = panoply_cna_setup.outputs,
      type = ome_type,
      subType = data_sub_type,
      fdr_cna_corr = fdr_for_cna_corr,
      params = add_parameters
  }

  call cna_corr_report_wdl.pgdac_cna_correlation_report {
    input:
      tarball = panoply_cna_correlation.outputs,
      type = ome_type,
      label = job_identifier,
      fdr = fdr_for_reports,
      tmpDir = "tmp"
  }
  
  call assoc_wdl.panoply_association {
    input: 
      inputData = panoply_cna_correlation.outputs, 
      groupsFile = association_groups,
      type = ome_type,
      subType = data_sub_type,
      params = add_parameters,
      standalone = standalone
  }

  call accum_wdl.panoply_accumulate as accumulate_assoc {
    input:
      input_tar = panoply_association.outputs,
      module = "association"
  } 

  Array[File] list_gct_assoc = accumulate_assoc.list_gct
  scatter (f in list_gct_assoc){
    call ssgsea_wdl.pgdac_ssgsea as ssgsea_assoc {
      input:
        input_ds = "${f}",
        gene_set_database = geneset_db,
        output_prefix = job_identifier,
        level = "gc"
    }
  }

  call cons_clust_wdl.pgdac_cons_clust {
    input:
      tarball = panoply_association.outputs,
      type = ome_type,
      groupsFile = cluster_enrichment_groups,
      subType = data_sub_type,
      params = add_parameters
  }

  call accum_wdl.panoply_accumulate as accumulate_clustering {
    input:
      input_tar = pgdac_cons_clust.outputs,
      module = "clustering"
  }

  Array[File] list_gct_clustering = accumulate_clustering.list_gct
  scatter (f in list_gct_clustering){
    call ssgsea_wdl.pgdac_ssgsea as ssgsea_clustering {
      input:
        input_ds = "${f}",
        gene_set_database = geneset_db,
        output_prefix = job_identifier,
        level = "gc"
    }
  }

  call download_wdl.panoply_download {
    input:
      cons_clust_tar = pgdac_cons_clust.outputs,
      ssgsea_ome_tar = ssgsea_ome.results,
      ssgsea_rna_tar = ssgsea_rna.results,
      analysisDir = job_identifier,
      ssgsea_assoc_tars = ssgsea_assoc.results,
      ssgsea_clust_tars = ssgsea_clustering.results,
      ptmsea = ptmsea_ome.results
  }

  output {
    File summary_and_ssgsea=panoply_download.summary
    File pgdac_full = panoply_download.full
    File rna_corr_report = pgdac_rna_protein_correlation_report.report
    File cna_corr_report = pgdac_cna_correlation_report.report
    File sample_qc_report = pgdac_sampleqc_report.report
  }

}
