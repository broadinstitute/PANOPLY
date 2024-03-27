version development
workflow panoply_spectronaut_libgen { 
  call spectronaut_libgen
}

task spectronaut_libgen {
  input {
    String license_key

    String experiment_name
    File? analysis_settings
    File? pulsar_settings
    # search databases -- upto 2 fasta files can be provided
    File fasta       
    File? fasta_1
    # spectral libraries -- upto 2 can be provided; if none specified, perform DirectDIA
    File? search_archive
    File? search_archive_1

    Directory files_folder
    File? file_of_files

    Int num_preemptions=0
    Int num_cpus=32
    Int ram_gb=512
    Int local_disk_gb=4000
  }

  Array[File] files = if defined(file_of_files) then read_lines(select_first([file_of_files])) else []
  Boolean directory_input = if defined(file_of_files) then false else true
  String raw_files = if directory_input then files_folder else sep(' -r ', files)
  String spec_lib = "${experiment_name}-speclib.kit"
  String search_archive_out = "${experiment_name}-archive.psar"
  
  
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
    /usr/bin/spectronaut -lg -se Pulsar ${"-es " + analysis_settings} ${"-rs " + pulsar_settings} \
        -n ${experiment_name} -a "$out_dir/${search_archive_out}" -k "$out_dir/${spec_lib}" \
        -fasta ${fasta} ${"-fasta " + fasta_1} ${"-sa " + search_archive} ${"-sa " + search_archive_1} \
        -d data -setTemp $sn_temp
    spectronaut -deactivate

    zip -r $out_zip $out_dir -x \*.zip

    # directory structure after completion of run ($working_dir = /root):
    #   /root/data/                                         input data
    #        /$out_dir/[timestamp]_${spec_lib}              spectral library (output)
    #        /$out_dir/${search_archive_out}                search archive output
    #        /$out_zip                                      output zip file

    cp $out_dir/*_${spec_lib} /$cromwell_root/${spec_lib}
    cp $out_dir/${search_archive_out} /$cromwell_root/${search_archive_out}
    mv $out_zip /$cromwell_root/
    
  }

  output {
    File spectronaut_output="spectronaut_output.zip"
    File spectral_library="${spec_lib}"
    File search_archive="${search_archive_out}"
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


