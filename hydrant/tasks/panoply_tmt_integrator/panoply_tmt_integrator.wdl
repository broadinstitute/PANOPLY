task panoply_tmt_integrator_task {
  File philosopher_output
  Float? ram_gb
  Int? local_disk_gb
  Int? num_preemptions
  Int? num_cpus

  File params
  File database
   
  command <<<
    set -euo pipefail
    
    unzip ${philosopher_output}
    
    awk '{ gsub( "protein_database:", "protein_database: ${database}")}1' \
      ${params} > tmt_int.yml
    
    collect_psm.sh ${philosopher_output}
    java -Xmx${if defined(ram_gb) then (ram_gb) else '52'}G \
      -jar /TMTIntegrator_v1.1.0.jar tmt_int.yml psms/*_psm.tsv
    zip results.zip /TMTIntegrator_v1.1.0/results/*
    
  >>>

  output {
    File results="results.zip"
  }

  runtime {
    docker      : "broadcptac/panoply_tmt_integrator:1"
    memory      : "${if defined(ram_gb) then ram_gb else '60'}GB"
    disks       : "local-disk ${if defined(local_disk_gb) then local_disk_gb else '100'} HDD"
    preemptible : "${if defined(num_preemptions) then num_preemptions else '0'}"
    cpu         : "${if defined(num_cpus) then num_cpus else '16'}"
  }

  meta {
    author: "Felipe da Veiga Leprevost"
    email : "felipe@leprevost.com.br"
  }
}

workflow panoply_tmt_integrator {
  Boolean package
  String null_file="gs://broad-institute-gdac/GDAC_FC_NULL"
  File database

  call panoply_tmt_integrator_task {
    input: database=database
  }
}

