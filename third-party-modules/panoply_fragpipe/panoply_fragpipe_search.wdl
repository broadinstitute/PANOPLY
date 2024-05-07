version development
workflow panoply_fragpipe_search { 
  input {
    File fragpipe_workflow
    File database
    String files_folder

    File? file_of_files
    String raw_file_type="DDA"
    File? fragpipe_manifest
    Boolean try_one_file=false
    
    Int? num_preemptions
    Int? num_cpus
    Int? ram_gb
    Int? local_disk_gb
  }

  call fragpipe {
    input: 
      fragpipe_workflow=fragpipe_workflow,
      database=database,
      files_folder=files_folder,

      file_of_files=file_of_files,
      fragpipe_manifest=fragpipe_manifest,
      raw_file_type=raw_file_type,
      try_one_file=try_one_file,

      num_preemptions=num_preemptions,
      num_cpus=num_cpus,
      ram_gb=ram_gb,
      local_disk_gb=local_disk_gb
  }
}


task fragpipe {
  input {
    File fragpipe_workflow
    File database
    
    String files_folder
    File? file_of_files
    String raw_file_type

    File? fragpipe_manifest
    Boolean try_one_file

    Int num_preemptions=0
    Int num_cpus=32
    Int ram_gb=128
    Int local_disk_gb=2000
  }
  Array[File] files = if defined(file_of_files) then read_lines(select_first([file_of_files])) else []

  command {
    . /etc/profile
    set -euo pipefail

    projdir="fragpipe"
    proc_data_zip="fragpipe_processed_data.zip"
    out_zip="fragpipe_output.zip"
    cromwell_root=$(pwd)                           # use cromwell_root fs for working dir
    working_dir=$(mktemp -d working_dir_XXXXXX)    # use wd in the /cromwell_root file system
    cd $working_dir
    
    mkdir -p $projdir/out
    chmod -R 777 $projdir

    cp -s ~{fragpipe_workflow} $projdir/
    cp -s ~{database} $projdir/
    frag_workflow=$(basename ${fragpipe_workflow})

    cd $projdir
    if [ -z ~{file_of_files} ]
    then
      tmp_dir=$(mktemp -d data_XXXXXX)
      mv ~{files_folder}/* $tmp_dir 
      mv $tmp_dir data
    else
      mkdir data
      cp ~{sep(' ', files)} data 
    fi

    echo "database.db-path=~{basename(database)}" >> $frag_workflow

    if [ "${try_one_file}" = true ]
    then
      file_to_keep=$(ls -1 "data" | head -1)
      cd data && ls | grep -v $file_to_keep | xargs rm -r && cd ..
    fi

    if [ -z ~{fragpipe_manifest} ]
    then
      python /usr/local/bin/get_fp_manifest.py /cromwell_root/$working_dir/$projdir/"data" "data" ~{raw_file_type} 
      frag_manifest="generated.fp-manifest"
    else
      cp -s ~{fragpipe_manifest} .
      frag_manifest=$(basename ${fragpipe_manifest})
      sed -i -e "s/\/path\//\/cromwell_root\/$working_dir\/$projdir\/data\//g" $frag_manifest  
    fi
    
    #headless version 
    /fragpipe/bin/fragpipe --headless --workflow $frag_workflow --manifest $frag_manifest --workdir "out" --config-msfragger /MSFragger-3.8/MSFragger-3.8.jar --config-philosopher /usr/local/bin/philosopher --config-ionquant /IonQuant-1.9.8/IonQuant-1.9.8.jar --config-python /opt/conda/envs/fragpipe/bin/python

    cd ..
    zip -r $out_zip $projdir/out -x \*.zip
    zip -r $proc_data_zip $projdir/data -x \*.zip

    mv $out_zip /$cromwell_root/
    mv $proc_data_zip /$cromwell_root/
  }

  output {
    File fragpipe_output="fragpipe_output.zip"
    File fragpipe_processed_data="fragpipe_processed_data.zip"
  }

  runtime {
    docker: "broadcptacdev/panoply_fragpipe:latest"
    cpuPlatform : "AMD Rome"
    memory: "${ram_gb}GB"
    bootDiskSizeGb: 512
    disks : "local-disk ${local_disk_gb} HDD"
    preemptible : num_preemptions
    cpu: num_cpus
  }

  meta {
    author: "Khoi Pham Munchic"
    email : "kpham@broadinstitute.org"
  }
}
