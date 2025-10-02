#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
task panoply_unified_assemble_results {
  ## main pipeline
  Array[File?] main_full
  Array[File?] main_summary
  Array[File?] norm_report
  Array[File?] rna_corr_report
  Array[File?] cna_corr_report
  Array[File?] ssgsea_ome_report
  Array[File?] ptmsea_ome_report
  Array[File?] omicsev_report
  Array[File?] cosmo_report
  Array[File?] sampleqc_report
  Array[File?] assoc_report
  Array[File?] blacksheep_report
  Array[File?] cmap_output
  Array[File?] cmap_ssgsea_output

  ## rna main
  File? ssgsea_rna_report
  File? rna_blacksheep_report
  File? rna_assoc_report
  File? immune_tar
  File? immune_report

  ## unified pipeline
  Array[File?]? clumpsptm_results
  File? clumpsptm_report
  Array[File?]? metaboanalyst_results
  Array[File?]? metaboanalyst_reports
  File? nmf_results
  File? nmf_reports

  String output_results_zip = "all_results.zip"
  String output_reports_zip = "all_reports.zip"

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions

  command {
    set -euo pipefail

    ### Setup RESULTS and REPORTS directory structure
    mkdir results
    mkdir results/proteogenomics_analysis results/nmf results/rna_analysis results/clumpsptm results/metaboanalyst
    mkdir reports
    mkdir reports/proteogenomics_analysis reports/nmf reports/rna_analysis reports/clumpsptm reports/metaboanalyst

    ### Dump results files into the given folders
    # MAIN 
    if [ ${sep='' main_full} != '' ]; then
      mv ${sep=' ' main_full} results/proteogenomics_analysis
      for filename in results/proteogenomics_analysis/*.tar;do tar -C results/proteogenomics_analysis -xvf $filename;rm $filename;done
    fi

    if [ ${sep='' main_summary} != '' ]; then
      mkdir results/proteogenomics_analysis/summary_files
      mv ${sep=' ' main_summary} results/proteogenomics_analysis/summary_files
      for filename in results/proteogenomics_analysis/summary_files/*.tar;do tar -C results/proteogenomics_analysis/summary_files -xvf $filename;rm $filename;done
    fi
    
    if [ ${sep='' cmap_output} != '' ]; then
      mkdir results/proteogenomics_analysis/proteome_cmap_analysis
      mv ${sep=' ' cmap_output} results/proteogenomics_analysis/proteome_cmap_analysis
      for filename in results/proteogenomics_analysis/proteome_cmap_analysis/*.tar;do tar -C results/proteogenomics_analysis/proteome_cmap_analysis -xvf $filename;rm $filename;done
    fi
    
    if [ ${sep='' cmap_ssgsea_output} != '' ]; then
      mv ${sep=' ' cmap_ssgsea_output} results/proteogenomics_analysis/proteome_cmap_analysis
      for filename in results/proteogenomics_analysis/proteome_cmap_analysis/*.tar;do tar -C results/proteogenomics_analysis/proteome_cmap_analysis -xvf $filename;rm $filename;done
    fi

    # MAIN REPORTS
    mkdir results/proteogenomics_analysis/all_html_reports
    if [ ${sep='' norm_report} != '' ]; then
      cp ${sep=' ' norm_report} results/proteogenomics_analysis/all_html_reports
      mv ${sep=' ' norm_report} reports/proteogenomics_analysis
    fi

    if [ ${sep='' rna_corr_report} != '' ]; then
      cp ${sep=' ' rna_corr_report} results/proteogenomics_analysis/all_html_reports
      mv ${sep=' ' rna_corr_report} reports/proteogenomics_analysis
    fi

    if [ ${sep='' cna_corr_report} != '' ]; then
      cp ${sep=' ' cna_corr_report} results/proteogenomics_analysis/all_html_reports
      mv ${sep=' ' cna_corr_report} reports/proteogenomics_analysis
    fi

    if [ ${sep='' ssgsea_ome_report} != '' ]; then
      cp ${sep=' ' ssgsea_ome_report} results/proteogenomics_analysis/all_html_reports
      mv ${sep=' ' ssgsea_ome_report} reports/proteogenomics_analysis
    fi

    if [ ${sep='' ptmsea_ome_report} != '' ]; then
      cp ${sep=' ' ptmsea_ome_report} results/proteogenomics_analysis/all_html_reports
      mv ${sep=' ' ptmsea_ome_report} reports/proteogenomics_analysis
    fi

    if [ ${sep='' omicsev_report} != '' ]; then
      cp ${sep=' ' omicsev_report} results/proteogenomics_analysis/all_html_reports
      mv ${sep=' ' omicsev_report} reports/proteogenomics_analysis
    fi
    
    if [ ${sep='' cosmo_report} != '' ]; then
      cp ${sep=' ' cosmo_report} results/proteogenomics_analysis/all_html_reports
      mv ${sep=' ' cosmo_report} reports/proteogenomics_analysis
    fi

    if [ ${sep='' sampleqc_report} != '' ]; then
      cp ${sep=' ' sampleqc_report} results/proteogenomics_analysis/all_html_reports
      mv ${sep=' ' sampleqc_report} reports/proteogenomics_analysis
    fi

    if [ ${sep='' assoc_report} != '' ]; then
      cp ${sep=' ' assoc_report} results/proteogenomics_analysis/all_html_reports
      mv ${sep=' ' assoc_report} reports/proteogenomics_analysis
    fi
    
    if [ ${sep='' blacksheep_report} != '' ]; then
      cp ${sep=' ' blacksheep_report} results/proteogenomics_analysis/all_html_reports
      mv ${sep=' ' blacksheep_report} reports/proteogenomics_analysis
    fi


    # RNA RESULTS
    mkdir results/rna_analysis/all_html_reports # make folder for all reports
    if [ ${immune_tar} != '' ]; then
      mkdir results/rna_analysis/immune_analysis
      mv ${immune_tar} results/rna_analysis/immune_analysis
      for filename in results/rna_analysis/immune_analysis/*.tar;do tar -C results/rna_analysis/immune_analysis -xvf $filename;rm $filename;done
    fi
    if [ ${immune_report} != '' ]; then
      cp ${immune_report} results/rna_analysis/all_html_reports
      mv ${immune_report} reports/rna_analysis
    fi
    if [ ${ssgsea_rna_report} != '' ]; then
      cp ${ssgsea_rna_report} results/rna_analysis/all_html_reports
      mv ${ssgsea_rna_report} reports/rna_analysis
    fi
    if [ ${rna_blacksheep_report} != '' ]; then
      cp ${rna_blacksheep_report} results/rna_analysis/all_html_reports
      mv ${rna_blacksheep_report} reports/rna_analysis
    fi
    if [ ${rna_assoc_report} != '' ]; then
      cp ${rna_assoc_report} results/rna_analysis/all_html_reports
      mv ${rna_assoc_report} reports/rna_analysis
    fi

    # UNIFIED RESULTS

    # NMF Results
    if [ ${nmf_results} != '' ]; then
      tar -C results/nmf -xvf ${nmf_results} --strip-components 1 # note: results tar already contains reports
    fi
    if [ ${nmf_reports} != '' ]; then
      tar -C reports/nmf -xvf ${nmf_reports} --strip-components 1
    fi

    # ClumpsPTM
    if [ ${sep='' clumpsptm_results} != '' ]; then
      mv ${sep=' ' clumpsptm_results} results/clumpsptm
      for filename in results/clumpsptm/*.tar;do 
        foldername=$(basename "$filename" .tar)
        mkdir -p "results/clumpsptm/$foldername"
        tar -C "results/clumpsptm/$foldername" -xvf "$filename"
        rm "$filename"
      done
    fi
    if [ ${clumpsptm_report} != '' ]; then
      cp ${clumpsptm_report} results/clumpsptm/all_html_reports
      mv ${clumpsptm_report} reports/clumpsptm
    fi

    # MetaboAnalyst
    if [ ${sep='' metaboanalyst_results} != '' ]; then
      mv ${sep=' ' metaboanalyst_results} results/metaboanalyst
      for filename in results/metaboanalyst/*.tar.gz;do 
        foldername=$(basename "$filename" .tar.gz)
        mkdir -p "results/metaboanalyst/$foldername"
        tar -C "results/metaboanalyst/$foldername" -xvf "$filename"
        rm "$filename"
      done
    fi
    if [ ${sep='' metaboanalyst_reports} != '' ]; then
      cp ${sep=' ' metaboanalyst_reports} results/metaboanalyst
      mv ${sep=' ' metaboanalyst_reports} reports/metaboanalyst
    fi

    ### Zip up final directories
    zip ${output_results_zip} -r results
    zip ${output_reports_zip} -r reports

  }

  output {
    File all_results = "${output_results_zip}"
    File all_reports = "${output_reports_zip}"

  }
  
  runtime {
    docker : "broadcptacdev/panoply_common:latest"
    memory : select_first ([memory, 16]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 64]) + " SSD"
    cpu : select_first ([num_threads, 1]) + ""
    preemptible : select_first ([num_preemptions, 0])
  }

  meta {
    author : "Myranda Maynard"
    email : "proteogenomics@broadinstitute.org"
  }
}

workflow panoply_unified_assemble_results_workflow {
  call panoply_unified_assemble_results
}