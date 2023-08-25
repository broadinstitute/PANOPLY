#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_normalize_filter_workflow/versions/14/plain-WDL/descriptor" as norm_filt_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_main_NG/versions/9/plain-WDL/descriptor" as main_NG_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_blacksheep_workflow/versions/13/plain-WDL/descriptor" as blacksheep_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_so_nmf_workflow/versions/33/plain-WDL/descriptor" as so_nmf_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_mo_nmf_gct/versions/13/plain-WDL/descriptor" as mo_nmf_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_so_nmf_sankey_workflow/versions/12/plain-WDL/descriptor" as so_nmf_sankey_wdl


workflow panoply_unified_NG_workflow {
  File? prote_ome
  File? phospho_ome
  File? acetyl_ome
  File? ubiquityl_ome
  File? sample_annotation
  File yaml
  String job_id
  String run_cmap
  String run_nmf #'true' or 'false'
  String? run_ptmsea

  # Normalize specific optional params:
  String? normalizeProteomics # "true" or "false"
  String? filterProteomics # "true" or "false"

  # Organize omics data into pairs

  Array[Pair[String?, File?]] ome_pairs =
    [ ("proteome", prote_ome),
      ("phosphoproteome", phospho_ome),
      ("acetylome", acetyl_ome),
      ("ubiquitylome", ubiquityl_ome) ]

  ### NORMALIZE:
  ### Normalize the data first so downstream modules (NMF etc) can run in parallel to main:
  scatter (pair in ome_pairs) {
    if ("${pair.right}" != "") {
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
  }
  
  # Zip Normalization / Filter Output into a labelled Array 
  Array[Pair[String?,File?]] ome_pairs_norm_filt = zip(norm_filt.output_ome_type , norm_filt.filtered_data_table)


  ### MAIN:
  scatter (pair in ome_pairs_norm_filt) {
    if ("${pair.right}" != "") {
        call main_NG_wdl.panoply_main_NG as pome {
          input:
            ## include all required arguments from above
            input_pome="${pair.right}",
            ome_type="${pair.left}",
            job_identifier="${job_id}-${pair.left}",
            run_ptmsea="${run_ptmsea}",
            run_nmf = "false",
            sample_annotation="${sample_annotation}",
            yaml="${yaml}"
        }
    }
  }
  
  ### BLACKSHEEP:
  scatter (pair in ome_pairs_norm_filt) {
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
      prote_ome = norm_filt.filtered_data_table[0],
      phospho_ome = norm_filt.filtered_data_table[1],
      acetyl_ome = norm_filt.filtered_data_table[2],
      ubiquityl_ome = norm_filt.filtered_data_table[3],
      run_sankey = "false" # run sankey_workflow separately
  }

  ### Multi-omics NMF:
  if ( run_nmf == "true" ){
    call mo_nmf_wdl.panoply_mo_nmf_gct_workflow as mo_nmf {
      input:
        yaml_file = yaml,
        label = job_id,
        omes = norm_filt.filtered_data_table
    }
  }
  
  ### NMF Sankey Diagrams (SO and MO nmf)
  call so_nmf_sankey_wdl.panoply_so_nmf_sankey_workflow as all_nmf_sankey {
  input:
    so_nmf_tar = so_nmf.nmf_results,
    mo_nmf_tar = mo_nmf.nmf_clust, #will exist if mo_nmf was run
    label = job_id
  }
  
  output {
    Array[File?] main_association_contrasts = pome.association_contrasts
    Array[File?] main_association_reports = pome.association_report
    
    Array[File?] main_ssgsea_results = pome.ssgsea_results
    Array[File?] main_ssgsea_reports = pome.ssgsea_report
    
    Array[File?]? main_ptmsea_results = pome.ptmsea_results
    Array[File?]? main_ptmsea_reports = pome.ptmsea_report

    Array[File?] blacksheep_tar = outlier.blacksheep_tar
    Array[File?] blacksheep_report = outlier.blacksheep_report

    File? mo_nmf_tar = mo_nmf.nmf_clust
    File? mo_nmf_report = mo_nmf.nmf_clust_report
    File? mo_nmf_ssgsea_report = mo_nmf.nmf_ssgsea_report

    File? so_nmf_results = so_nmf.nmf_results
    File? so_nmf_reports = so_nmf.nmf_reports

    File? so_nmf_sankey_results = all_nmf_sankey.sankey_tar
    File? so_nmf_sankey_report = all_nmf_sankey.sankey_report
  }
 }