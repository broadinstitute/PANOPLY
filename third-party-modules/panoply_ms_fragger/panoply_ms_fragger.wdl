task panoply_ms_fragger {
  Boolean package
  String package_name
  String null_file
  String package_archive="${package_name}.zip"
  Float? ram_gb
  Int? local_disk_gb
  Int? num_preemptions
  Int? num_cpus

  File params
  File database
  Array[File] mzML_files
  Float heap_fraction=0.875

  command <<<
    set -euo pipefail

    # set user-defined parameters
    awk '$1 == "database_name" {$3 = "${database}"}
         $1 == "num_threads"   {$3 = 0}
         {print}' ${params} > fragger.params

    # Run MSFragger
    java -Xmx${if defined(ram_gb) then floor(ram_gb * heap_fraction) else '52'}G \
      -jar /MSFragger-2.1.jar \
      fragger.params \
      ${sep=' ' mzML_files}

     zip -r pepXML . -i \*.pepXML

     if ${package}; then
       package.sh -x broad-institute-gdac/\* \
                  -x \*.mzML \
                  -x \*.pepXML \
                  -x pepXML.zip \
                  ${package_archive}
        fi
  >>>

  output {
    File printfile="${params}"
    File philosopher_pipeline_pkg="${if package then package_archive else null_file}"
    File pepXML_archive="pepXML.zip"
  }

  runtime {
    docker      : "broadcptac/panoply_ms_fragger:1"
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

workflow panoply_ms_fragger_workflow
{
  call panoply_ms_fragger
  {
    input: 
      package_name="test-ms-fragger",
      null_file="gs://broad-institute-gdac/GDAC_FC_NULL"
  }
}
