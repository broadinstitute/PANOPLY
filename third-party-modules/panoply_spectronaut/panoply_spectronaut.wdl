version development
workflow panoply_spectronaut { 
  call spectronaut
}

task spectronaut {
  input {
    String license_key="b233d39d-b28c-47b9-8e41-e16a9ba923ad"

    String experiment_name
    File? analysis_settings
    File? condition_setup
    # search databases -- upto 2 fasta files can be provided
    File fasta       
    File? fasta_1
    # spectral libraries -- upto 2 can be provided; if none specified, perform DirectDIA
    File? spectral_library
    File? spectral_library_1
    File? report_schema
    File? json_settings

    Directory files_folder
    File? file_of_files

    Int num_preemptions=0
    Int num_cpus=32
    Int ram_gb=512
    Int local_disk_gb=4000
  }

  Array[File] files = if defined(file_of_files) then read_lines(select_first([file_of_files])) else []
  Boolean directory_input = if defined(file_of_files) then false else true
  Boolean direct_DIA = if defined(spectral_library) then false else true
  String raw_files = if directory_input then files_folder else sep(' -r ', files)
  
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
    spectronaut -activate ${license_key}
    /usr/bin/spectronaut ${if direct_DIA then "-direct" else ""} ${"-s " + analysis_settings} \
        ${"-con " + condition_setup} -n ${experiment_name} -o $out_dir \
        -fasta ${fasta} ${"-fasta " + fasta_1} ${"-a " + spectral_library} ${"-a " + spectral_library_1} \
        ${"-rs " + report_schema} ${"-j " + json_settings} -d data -setTemp $sn_temp
    spectronaut -deactivate

    zip -r $out_zip $out_dir -x \*.zip

    # directory structure after completion of run ($working_dir = /root):
    #   /root/data/                                         input data
    #        /$out_dir/[timestamp]_${experiment_name}/      output
    #        /$out_zip                                      output zip file

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
    author: "D. R. Mani"
    email : "proteogenomics@broadinstitute.org"
  }
}


