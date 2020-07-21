task panoply_philosopher {
  Boolean package
  String null_file
  File package_archive
  String package_name=basename(package_archive)
  Float? ram_gb
  Int? local_disk_gb
  Int? num_preemptions
  Int? num_cpus

  File params
  File database
  Array[File] mzML_files
  File pepXML_archive
  File mapping_table
  String output_prefix

  command <<<
    . /etc/profile

    set -euo pipefail

    unzip ${pepXML_archive}
    awk '{ gsub( "protein_database:", "protein_database: ${database}")}1' \
      ${params} > philosopher.yml

    sed -r -i -e 's/Proteome/W/' \
              -e 's/Phosphoproteome/P/' \
              -e 's/\b([[:digit:]]+)(CPTAC.*)POOL/\1\2pool\1/' \
              -e 's/\b131\b/&N/g' ${mapping_table}

    rootdir=$(dirname `dirname ${select_first( mzML_files )}`)

    # Symlink project directory
    projdir=${output_prefix}_philosopher
    ln -s $rootdir $projdir

    # Go to project directory for generating annotation files and running
    # Philosopher
    mv philosopher.yml $projdir/philosopher.yml
    cd $projdir

    # Generate annotation files from modified mapping file
    awk '
         NR == 1 {
            for (i=3; i<=NF; i++) {
              if ($i ~ /^[0-9]+/) {
                channels[i] = $i
              }
            }
          }
          NR > 1 {
            for (i=3; i<=NF; i++) {
              if (i in channels) {
                print channels[i], $i > $2"/annotation.txt"
              }
            }
          }' ${mapping_table}

    # Run Philosopher
    philosopher pipeline --config philosopher.yml */

    # Archive results
    cd -
    zip -r $projdir $projdir -x \*.mzML

    if ${package}; then
      mv ${package_archive} .
      package.sh -x broad-institute-gdac/\* \
                 -x \*.mzML \
                 -x $projdir.zip \
                 ${package_name}
    fi
 
  >>>

  output {
    File philosopher_pipeline_pkg="${if package then package_name else null_file}"
    File philosopher_results="${output_prefix}_philosopher.zip"
  }

  runtime {
    docker : "broadcptac/panoply_ms_philosopher:1"
    memory: "${if defined(ram_gb) then ram_gb else '26'}GB"
    disks : "local-disk ${if defined(local_disk_gb) then local_disk_gb else '120'} HDD"
    preemptible : "${if defined(num_preemptions) then num_preemptions else '0'}"
    cpu: "${if defined(num_cpus) then num_cpus else '4'}"
  }

  meta{
    author: "Felipe da Veiga Leprevost"
    email : "felipe@leprevost.com.br"
  }
}

workflow panoply_ms_philosopher {
  Boolean package
  String null_file="gs://broad-institute-gdac/GDAC_FC_NULL"
  String package_name="philosopher_pipeline"
  Array[File] mzML_files
  File database
  File pepXML_archive
  File package_archive

  call panoply_philosopher {
    input: package=package,
           null_file=null_file,
           package_archive=package_archive,
           mzML_files=mzML_files,
           pepXML_archive=pepXML_archive,
           database=database
  }
}
