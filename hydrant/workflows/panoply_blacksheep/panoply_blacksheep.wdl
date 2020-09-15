task panoply_blacksheep {
  Int? memory
  Int? disk_space
  Int? num_threads
  
  File input_gct
  File master_yaml
  String output_prefix
  
  String? apply_filtering
  File? identifiers_file
  File? groups_file
  Float? fraction_samples_cutoff
  Float? fdr_value
  
  command {
    set -euo pipefail
    
    /usr/bin/Rscript /prot/proteomics/Projects/PGDAC/src/parameter_manager.r \
    --module blacksheep \
    --master_yaml ${master_yaml} \
    ${"--blacksheep_apply_filtering " + apply_filtering} \
    ${"--blacksheep_identifiers_file " + identifiers_file} \
    ${"--blacksheep_groups_file " + groups_file} \
    ${"--blacksheep_fraction_samples_cutoff " + fraction_samples_cutoff} \
    ${"--blacksheep_fdr_value " + fdr_value}
    
    /usr/bin/Rscript /prot/proteomics/Projects/PGDAC/src/blacksheep_rcode.R "${input_gct}" "final_output_params.yaml"
    
    if [ "${groups_file}" != "" ]; then
    cp ${groups_file} "blacksheep"
    fi
    
    tar -czvf "${output_prefix}_blacksheep.tar" blacksheep final_output_params.yaml
  }
  
  output {
    File tar_out = "${output_prefix}_blacksheep.tar"
  }
  
  runtime {
    docker : "broadcptacdev/panoply_blacksheep:latest"
    memory : select_first ([memory, 10]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 20]) + " SSD"
    cpu : select_first ([num_threads, 1]) + ""
  }
  
  meta {
    author : "Karen Christianson"
    email : "karen@broadinstitute.org"
  }
}

task panoply_blacksheep_report {
  File input_tar
  String output_prefix
  
  Int? memory
  Int? disk_space
  Int? num_threads
  
  command {
    set -euo pipefail
    
    /usr/bin/Rscript /home/pgdac/src/rmd_blacksheep.R "${input_tar}" "${output_prefix}"
  }
  
  output {
    File report_out = "${output_prefix}_blacksheep_rmd.html"
  }
  
  runtime {
    docker : "broadcptacdev/panoply_blacksheep_report:latest"
    memory : select_first ([memory, 32]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 64]) + " SSD"
    cpu : select_first ([num_threads, 1]) + ""
  }
  
  meta {
    author : "Karen Christianson"
    email : "karenchristianson@broadinstitute.org"
  }
}

workflow panoply_blacksheep_workflow {
  File input_gct
  File master_yaml
  String output_prefix
  File? groups_file
  
  call panoply_blacksheep {
    input:
      input_gct = input_gct,
      master_yaml = master_yaml,
      output_prefix = output_prefix,
      groups_file = groups_file
  }
  
  call panoply_blacksheep_report {
    input:
      input_tar = panoply_blacksheep.tar_out,
      output_prefix = output_prefix
  }
  
  output{
    File blacksheep_tar = panoply_blacksheep.tar_out
    File blacksheep_report = panoply_blacksheep_report.report_out
  }
  
}