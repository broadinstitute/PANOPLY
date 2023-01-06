version development
workflow panoply_fragpipe_search { 
  input {
    File fragpipe_workflow
    File database
    Directory files_folder

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
    
    Directory files_folder
    File? file_of_files
    String raw_file_type

    File? fragpipe_manifest
    Boolean try_one_file

    Int? num_preemptions
    Int? num_cpus
    Int? ram_gb
    Int? local_disk_gb
  }
  Array[File] files = if defined(file_of_files) then read_lines(select_first([file_of_files])) else []

  command {
    . /etc/profile
    set -euo pipefail

    projdir="fragpipe"
    proc_data_zip="fragpipe_processed_data.zip"
    out_zip="fragpipe_output.zip"
    cromwell_root=$(pwd)
    cd ~
    mkdir -p $projdir/out
    chmod -R 777 $projdir

    cp -s ~{fragpipe_workflow} $projdir/
    cp -s ~{database} $projdir/
    frag_workflow=$(basename ${fragpipe_workflow})

    cd $projdir
    if [ -z ~{file_of_files} ]
    then
      cp -a ~{files_folder} .
      mv $(basename ~{files_folder}) data
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
      python /usr/local/bin/get_fp_manifest.py /root/$projdir/"data" "data" ~{raw_file_type}
      frag_manifest="generated.fp-manifest"
    else
      cp -s ~{fragpipe_manifest} .
      frag_manifest=$(basename ${fragpipe_manifest})
    fi
    
    /fragpipe/bin/fragpipe --headless --workflow $frag_workflow --manifest $frag_manifest --workdir "out" --config-msfragger /MSFragger-3.6/MSFragger-3.6.jar --config-philosopher /usr/local/bin/philosopher --config-ionquant /IonQuant-1.8.9/IonQuant-1.8.9.jar --config-python /opt/conda/envs/fragpipe/bin/python

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
    memory: "${if defined(ram_gb) then ram_gb else '128'}GB"
    bootDiskSizeGb: 512
    disks : "local-disk ${if defined(local_disk_gb) then local_disk_gb else '2000'} HDD"
    preemptible : "${if defined(num_preemptions) then num_preemptions else '0'}"
    cpu: "${if defined(num_cpus) then num_cpus else '32'}"
  }

  meta {
    author: "Khoi Pham Munchic"
    email : "kpham@broadinstitute.org"
  }
}
