version development
workflow panoply_spectronaut_SNE { 
  call spectronaut_SNE
}

task spectronaut_SNE {
  input {

    String experiment_name
    File fasta
    Directory files_folder
    Boolean sne_out=false

    # Specify an additional FASTA
    File? fasta_1
    File? settings_schema

    # Specify report schema outputs
    File? report_schema
    File? report_schema_1
    File? report_schema_2

   
    File? file_of_files

    Int num_preemptions=0
    Int num_cpus=32
    Int ram_gb=512
    Int local_disk_gb=4000
  }

  Array[File] files = if defined(file_of_files) then read_lines(select_first([file_of_files])) else []
  Boolean directory_input = if defined(file_of_files) then false else true
  String sne_files = if directory_input then files_folder else sep(' -sne ', files)
  
  command {
    set -euo pipefail

    out_zip="spectronaut_output.zip"
    out_dir="spectronaut/out"
    cromwell_root=$(pwd)                           # use cromwell_root fs for both wd and temp dir
    sn_temp=$(mktemp -d sn_temp_XXXXXX)            # temp dir for Spectronaut -- else runs out of space on root fs 
    working_dir=$(mktemp -d working_dir_XXXXXX)    # use wd in the /cromwell_root file system
    cd $working_dir

    mkdir -p $out_dir
    if [[ "${directory_input}" = "true" ]]
    then
      tmp_dir=$(mktemp -d data_XXXXXX)      # in case the files_folder is named 'data'
      mv ${files_folder}/* $tmp_dir         # all under $cromwell_root -- no need to copy
      mv $tmp_dir data
    else
      mkdir data
      cp ${sep(' ', files)} data
    fi
    
    # run spectronaut
    if [[ "${sne_out}" = "true" ]]
    then
        /usr/bin/spectronaut manageSNE --merge -o $out_dir -d data -n ${experiment_name} \
          ${"-rs " + report_schema} ${"-rs " + report_schema_1} ${"-rs " + report_schema_2} -setTemp $sn_temp
    else
        /usr/bin/spectronaut combine -o $out_dir -d data \
          -fasta ${fasta} ${"-fasta " + fasta_1} ${"-s " + settings_schema} ${"-rs " + report_schema} \
          ${"-rs " + report_schema_1} ${"-rs " + report_schema_2} -n ${experiment_name} -setTemp $sn_temp
    fi

    zip -r $out_zip $out_dir -x \*.zip 

    mv $out_zip /$cromwell_root/
    
  }

  output {
    File spectronaut_output="spectronaut_output.zip"
  }

  runtime {
    docker: "broadcptacdev/panoply_spectronaut:latest"
    cpuPlatform: "AMD Rome"
    memory: "${ram_gb}GB"   # 896GB max for AMD Rome
    bootDiskSizeGb: 512
    disks : "local-disk ${local_disk_gb} HDD"
    preemptible : num_preemptions
    cpu: num_cpus
  }

  meta {
    author: "Simone Gohsman"
    email : "proteogenomics@broadinstitute.org"
  }
}