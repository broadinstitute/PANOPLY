#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
version 1.0

task panoply_unified_assemble_results {
  input {
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
  }

  command {
    set -euo pipefail

    ### Setup RESULTS and REPORTS directory structure
    mkdir results
    mkdir results/proteogenomics_analysis results/nmf results/rna_analysis results/clumpsptm results/metaboanalyst
    mkdir reports
    mkdir reports/proteogenomics_analysis reports/nmf reports/rna_analysis reports/clumpsptm reports/metaboanalyst

    ### Dump results files into the given folders
    # MAIN 
    if [ ${sep='' select_all(main_full)} != '' ]; then
      mv ${sep=' ' select_all(main_full)} results/proteogenomics_analysis
      for filename in results/proteogenomics_analysis/*.tar;do tar -C results/proteogenomics_analysis -xvf $filename;rm $filename;done
    fi

    if [ ${sep='' select_all(main_summary)} != '' ]; then
      mkdir results/proteogenomics_analysis/summary_files
      mv ${sep=' ' select_all(main_summary)} results/proteogenomics_analysis/summary_files
      for filename in results/proteogenomics_analysis/summary_files/*.tar;do tar -C results/proteogenomics_analysis/summary_files -xvf $filename;rm $filename;done
    fi
    
    if [ ${sep='' select_all(cmap_output)} != '' ]; then
      mkdir results/proteogenomics_analysis/proteome_cmap_analysis
      mv ${sep=' ' select_all(cmap_output)} results/proteogenomics_analysis/proteome_cmap_analysis
      for filename in results/proteogenomics_analysis/proteome_cmap_analysis/*.tar;do tar -C results/proteogenomics_analysis/proteome_cmap_analysis -xvf $filename;rm $filename;done
    fi
    
    if [ ${sep='' select_all(cmap_ssgsea_output)} != '' ]; then
      mv ${sep=' ' select_all(cmap_ssgsea_output)} results/proteogenomics_analysis/proteome_cmap_analysis
      for filename in results/proteogenomics_analysis/proteome_cmap_analysis/*.tar;do tar -C results/proteogenomics_analysis/proteome_cmap_analysis -xvf $filename;rm $filename;done
    fi

    # MAIN REPORTS
    mkdir results/proteogenomics_analysis/all_html_reports
    if [ ${sep='' select_all(norm_report)} != '' ]; then
      cp ${sep=' ' select_all(norm_report)} results/proteogenomics_analysis/all_html_reports
      mv ${sep=' ' select_all(norm_report)} reports/proteogenomics_analysis
    fi

    if [ ${sep='' select_all(rna_corr_report)} != '' ]; then
      cp ${sep=' ' select_all(rna_corr_report)} results/proteogenomics_analysis/all_html_reports
      mv ${sep=' ' select_all(rna_corr_report)} reports/proteogenomics_analysis
    fi

    if [ ${sep='' select_all(cna_corr_report)} != '' ]; then
      cp ${sep=' ' select_all(cna_corr_report)} results/proteogenomics_analysis/all_html_reports
      mv ${sep=' ' select_all(cna_corr_report)} reports/proteogenomics_analysis
    fi

    if [ ${sep='' select_all(ssgsea_ome_report)} != '' ]; then
      cp ${sep=' ' select_all(ssgsea_ome_report)} results/proteogenomics_analysis/all_html_reports
      mv ${sep=' ' select_all(ssgsea_ome_report)} reports/proteogenomics_analysis
    fi

    if [ ${sep='' select_all(ptmsea_ome_report)} != '' ]; then
      cp ${sep=' ' select_all(ptmsea_ome_report)} results/proteogenomics_analysis/all_html_reports
      mv ${sep=' ' select_all(ptmsea_ome_report)} reports/proteogenomics_analysis
    fi

    if [ ${sep='' select_all(omicsev_report)} != '' ]; then
      cp ${sep=' ' select_all(omicsev_report)} results/proteogenomics_analysis/all_html_reports
      mv ${sep=' ' select_all(omicsev_report)} reports/proteogenomics_analysis
    fi
    
    if [ ${sep='' select_all(cosmo_report)} != '' ]; then
      cp ${sep=' ' select_all(cosmo_report)} results/proteogenomics_analysis/all_html_reports
      mv ${sep=' ' select_all(cosmo_report)} reports/proteogenomics_analysis
    fi

    if [ ${sep='' select_all(sampleqc_report)} != '' ]; then
      cp ${sep=' ' select_all(sampleqc_report)} results/proteogenomics_analysis/all_html_reports
      mv ${sep=' ' select_all(sampleqc_report)} reports/proteogenomics_analysis
    fi

    if [ ${sep='' select_all(assoc_report)} != '' ]; then
      cp ${sep=' ' select_all(assoc_report)} results/proteogenomics_analysis/all_html_reports
      mv ${sep=' ' select_all(assoc_report)} reports/proteogenomics_analysis
    fi
    
    if [ ${sep='' select_all(blacksheep_report)} != '' ]; then
      cp ${sep=' ' select_all(blacksheep_report)} results/proteogenomics_analysis/all_html_reports
      mv ${sep=' ' select_all(blacksheep_report)} reports/proteogenomics_analysis
    fi


    # RNA RESULTS
    mkdir results/rna_analysis/all_html_reports # make folder for all reports
    if [ "${defined(immune_tar)}" = "true" ]; then
      mkdir results/rna_analysis/immune_analysis
      mv ${immune_tar} results/rna_analysis/immune_analysis
      for filename in results/rna_analysis/immune_analysis/*.tar;do tar -C results/rna_analysis/immune_analysis -xvf $filename;rm $filename;done
    fi
    if [ "${defined(immune_report)}" = "true" ]; then
      cp ${immune_report} results/rna_analysis/all_html_reports
      mv ${immune_report} reports/rna_analysis
    fi
    if [ "${defined(ssgsea_rna_report)}" = "true" ]; then
      cp ${ssgsea_rna_report} results/rna_analysis/all_html_reports
      mv ${ssgsea_rna_report} reports/rna_analysis
    fi
    if [ "${defined(rna_blacksheep_report)}" = "true" ]; then
      cp ${rna_blacksheep_report} results/rna_analysis/all_html_reports
      mv ${rna_blacksheep_report} reports/rna_analysis
    fi
    if [ "${defined(rna_assoc_report)}" = "true" ]; then
      cp ${rna_assoc_report} results/rna_analysis/all_html_reports
      mv ${rna_assoc_report} reports/rna_analysis
    fi

    # UNIFIED RESULTS

    # NMF Results
    if [ "${defined(nmf_results)}" = "true" ]; then
      tar -C results/nmf -xvf ${nmf_results} --strip-components 1 # note: results tar already contains reports
    fi
    if [ "${defined(nmf_reports)}" = "true" ]; then
      tar -C reports/nmf -xvf ${nmf_reports} --strip-components 1
    fi

    # ClumpsPTM
    if [ "${defined(clumpsptm_results)}" = "true" ]; then
      mv ${if defined(clumpsptm_results) then sep(' ', select_all(select_first([clumpsptm_results]))) else ""} results/clumpsptm
      for filename in results/clumpsptm/*.tar;do 
        foldername=$(basename "$filename" .tar)
        mkdir -p "results/clumpsptm/$foldername"
        tar -C "results/clumpsptm/$foldername" -xvf "$filename"
        rm "$filename"
      done
    fi
    if [ "${defined(clumpsptm_report)}" = "true" ]; then
      mkdir -p results/clumpsptm/all_html_reports
      cp ${clumpsptm_report} results/clumpsptm/all_html_reports
      mv ${clumpsptm_report} reports/clumpsptm
    fi

    # MetaboAnalyst
    if [ "${defined(metaboanalyst_results)}" = "true" ]; then
      mv ${if defined(metaboanalyst_results) then sep(' ', select_all(select_first([metaboanalyst_results]))) else ""} results/metaboanalyst
      for filename in results/metaboanalyst/*.tar.gz;do 
        foldername=$(basename "$filename" .tar.gz)
        mkdir -p "results/metaboanalyst/$foldername"
        tar -C "results/metaboanalyst/$foldername" -xvf "$filename"
        rm "$filename"
      done
    fi
    if [ "${defined(metaboanalyst_reports)}" = "true" ]; then
      cp ${if defined(metaboanalyst_reports) then sep(' ', select_all(select_first([metaboanalyst_reports]))) else ""} results/metaboanalyst
      mv ${if defined(metaboanalyst_reports) then sep(' ', select_all(select_first([metaboanalyst_reports]))) else ""} reports/metaboanalyst
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