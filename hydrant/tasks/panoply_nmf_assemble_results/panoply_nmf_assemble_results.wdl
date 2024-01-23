#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
task panoply_nmf_assemble_results {
  File? mo_nmf_results
  File? mo_nmf_figures
  File? mo_nmf_report
  File? mo_nmf_ssgsea_tar
  File? mo_nmf_ssgsea_report

  Array[File?] so_nmf_results
  Array[File?] so_nmf_figures
  Array[File?] so_nmf_report
  Array[File?] so_nmf_ssgsea_tar
  Array[File?] so_nmf_ssgsea_report

  File? sankey_tar
  File? sankey_report

  String output_results_tar = "nmf_results.tar"
  String output_reports_tar = "nmf_reports.tar"

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions

  command {
    set -euo pipefail

    ### Setup RESULTS and REPORTS directory structure
    mkdir nmf_results
    mkdir nmf_reports



    ## Compile Multiomic Results
    if [ ${mo_nmf_results} != '' ]; then
      mkdir nmf_results/mo_nmf # make general mo_nmf directory

      ## untar panoply_nmf tar
      mkdir nmf_results/mo_nmf/results
      tar -C nmf_results/mo_nmf/results -zxvf $mo_nmf_results #untar file in appropriate directory

      ## untar panoply_nmf_postprocessing tar
      mkdir nmf_results/mo_nmf/figures
      tar -C nmf_results/mo_nmf/figures -zxvf $mo_nmf_figures #untar file in appropriate directory
      
      ## copy in report
      cp ${mo_nmf_report} nmf_results/mo_nmf/ #copy report to results
      cp ${mo_nmf_report} nmf_reports/ #copy report to reports
      
      ## copy in ssgsea results
      mkdir nmf_results/mo_nmf/mo_nmf_ssgsea
      tar -C nmf_results/mo_nmf/mo_nmf_ssgsea -zxvf $mo_nmf_ssgsea #untar file in appropriate directory

      ## copy in ssgsea report
      cp ${mo_nmf_ssgsea_report} nmf_results/mo_nmf/ #copy reports to results
      cp ${mo_nmf_ssgsea_report} nmf_reports #copy reports to reports

    fi



    ## Compile Single-Ome Results
    if [ ${so_nmf_results} != '' ]; then
      mkdir nmf_results/so_nmf # make general so_nmf directory

      ## untar panoply_nmf tar
      if [ ${sep='' so_nmf_results} != '' ]; then
        mkdir nmf_results/so_nmf/results
        for tar in ${sep=' ' so_nmf_results}
        do
          cp $tar nmf_results/so_nmf/results #copy tars into folder
          for filename in nmf_results/so_nmf/results/*.tar
              do mkdir nmf_results/so_nmf/results/$(basename "$filename" .tar)
              tar -C nmf_results/so_nmf/results/$(basename "$filename" .tar) -xvf $filename
              rm $filename #untar and remove
          done
        done
      fi


      ## untar panoply_nmf_postprocess tar
      if [ ${sep='' so_nmf_figures} != '' ]; then
        mkdir nmf_results/so_nmf/figures
        for tar in ${sep=' ' so_nmf_figures}
        do
          cp $tar nmf_results/so_nmf/figures #copy tars into folder
          for filename in nmf_results/so_nmf/figures/*.tar
              do mkdir nmf_results/so_nmf/figures/$(basename "$filename" .tar)
              tar -C nmf_results/so_nmf/figures/$(basename "$filename" .tar) -xvf $filename
              rm $filename #untar and remove
          done
        done
      fi
      
      # copy in reports
      mkdir nmf_results/so_nmf/so_nmf_report
      mkdir nmf_reports/so_nmf
      if [ ${sep='' so_nmf_report} != '' ]; then
        cp ${sep=' ' so_nmf_report} nmf_results/so_nmf/so_nmf_report #copy reports to results
        cp ${sep=' ' so_nmf_report} nmf_reports/so_nmf #copy reports to reports
      fi
      
      # copy in ssgsea results
      mkdir nmf_results/so_nmf/so_nmf_ssgsea
      if [ ${sep='' so_nmf_ssgsea_tar} != '' ]; then
        cp ${sep=' ' so_nmf_ssgsea_tar} nmf_results/so_nmf/so_nmf_ssgsea #copy tars in
        for filename in nmf_results/so_nmf/so_nmf_ssgsea/*.tar.gz
          do mkdir nmf_results/so_nmf/so_nmf_ssgsea/$(basename "$filename" .tar.gz) #make directory
      tar -C nmf_results/so_nmf/so_nmf_ssgsea/$(basename "$filename" .tar.gz) -zxvf $filename #untar file
          rm $filename #remove tar
        done #untar
      fi

      # copy in ssgsea reports
      mkdir nmf_results/so_nmf/so_nmf_ssgsea_report
      mkdir nmf_reports/so_nmf_ssgsea
      if [ ${sep='' so_nmf_ssgsea_report} != '' ]; then
        cp ${sep=' ' so_nmf_ssgsea_report} nmf_results/so_nmf/so_nmf_ssgsea_report #copy reports to results
        cp ${sep=' ' so_nmf_ssgsea_report} nmf_reports/so_nmf_ssgsea #copy reports to reports
      fi
    fi



    ## Compile Sankey Results
    if [ ${sankey_tar} != '' ]; then
      mkdir nmf_results/sankey # make general mo_nmf directory

      ## untar sankey_tar
      mkdir nmf_results/sankey
      tar -C nmf_results/sankey -zxvf $sankey_tar #untar file in appropriate directory

      ## copy in report
      cp ${sankey_report} nmf_results/sankey/ #copy report to results
      cp ${sankey_report} nmf_reports/ #copy report to reports
      
    fi



    ### Tar final directories
    tar -czf ${output_results_tar} nmf_results/
    tar -czf ${output_reports_tar} nmf_reports/
  }

  output {
    File nmf_results = "${output_results_tar}"
    File nmf_reports = "${output_reports_tar}"

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

workflow panoply_nmf_assemble_results_workflow {
  call panoply_nmf_assemble_results
}