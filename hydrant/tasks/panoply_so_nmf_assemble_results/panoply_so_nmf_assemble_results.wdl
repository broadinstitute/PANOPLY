#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
task panoply_so_nmf_assemble_results {
  Array[File?] so_nmf_tar
  Array[File?] so_nmf_report
  Array[File?] so_nmf_ssgsea_tar
  Array[File?] so_nmf_ssgsea_report

  String output_results_zip = "nmf_results.zip"
  String output_reports_zip = "nmf_reports.zip"

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions

  command {
    set -euo pipefail

    ### Setup RESULTS and REPORTS directory structure
    mkdir nmf_results
    mkdir nmf_reports
    # SO-NMF
    mkdir nmf_results/so-nmf
    if [ ${sep='' so_nmf_tar} != '' ]; then
      for tar in ${sep=' ' so_nmf_tar}
      do
        cp $tar nmf_results/so-nmf #copy tar in
        for filename in nmf_results/so-nmf/*.tar
            do mkdir nmf_results/so-nmf/$(basename "$filename" .tar)
            tar -C nmf_results/so-nmf/$(basename "$filename" .tar) -xvf $filename
            rm $filename #untar and remove
          done
      done
    fi
    
    mkdir nmf_results/so-nmf_report
    mkdir nmf_reports/so-nmf
    if [ ${sep='' so_nmf_report} != '' ]; then
      cp ${sep=' ' so_nmf_report} nmf_results/so-nmf_report #copy reports to results
      cp ${sep=' ' so_nmf_report} nmf_reports/so-nmf #copy reports to reports
    fi
    
    mkdir nmf_results/so-nmf_ssgsea
    if [ ${sep='' so_nmf_ssgsea_tar} != '' ]; then
      cp ${sep=' ' so_nmf_ssgsea_tar} nmf_results/so-nmf_ssgsea #copy tars in
      for filename in nmf_results/so-nmf_ssgsea/*.tar.gz
        do mkdir nmf_results/so-nmf_ssgsea/$(basename "$filename" .tar.gz) #make directory
    tar -C nmf_results/so-nmf_ssgsea/$(basename "$filename" .tar.gz) -zxvf $filename #untar file
        rm $filename #remove tar
      done #untar
    fi

    mkdir nmf_results/so-nmf_ssgsea_report
    mkdir nmf_reports/so-nmf_ssgsea
    if [ ${sep='' so_nmf_ssgsea_report} != '' ]; then
      cp ${sep=' ' so_nmf_ssgsea_report} nmf_results/so-nmf_ssgsea_report #copy reports to results
      cp ${sep=' ' so_nmf_ssgsea_report} nmf_reports/so-nmf_ssgsea #copy reports to reports
    fi

    ### Zip up final directories
    zip ${output_results_zip} -r nmf_results
    zip ${output_reports_zip} -r nmf_reports
  }

  output {
    File nmf_results = "${output_results_zip}"
    File nmf_reports = "${output_reports_zip}"

  }
  
  runtime {
    docker : "broadcptacdev/panoply_common:latest"
    memory : select_first ([memory, 16]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 32]) + " SSD"
    cpu : select_first ([num_threads, 1]) + ""
    preemptible : select_first ([num_preemptions, 0])
  }

  meta {
    author : "C. Williams"
    email : "proteogenomics@broadinstitute.org"
  }
}

workflow panoply_so_nmf_assemble_results_workflow {
  call panoply_so_nmf_assemble_results
}