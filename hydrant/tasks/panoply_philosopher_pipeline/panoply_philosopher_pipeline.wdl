task pgdac_msfragger 
{
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

  output 
  {
    File philosopher_pipeline_pkg="${if package then package_archive else null_file}"
    File pepXML_archive="pepXML.zip"
  }

  runtime 
  {
    docker      : "broadcptac/pgdac_ms_fragger:1"
    memory      : "${if defined(ram_gb) then ram_gb else '60'}GB"
    disks       : "local-disk ${if defined(local_disk_gb) then local_disk_gb else '100'} HDD"
    preemptible : "${if defined(num_preemptions) then num_preemptions else '0'}"
    cpu         : "${if defined(num_cpus) then num_cpus else '16'}"
  }

  meta 
  {
    author: "Felipe da Veiga Leprevost"
    email : "felipe@leprevost.com.br"
  }
}

task pgdac_philosopher 
{
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

  output 
  {
    File philosopher_pipeline_pkg="${if package then package_name else null_file}"
    File philosopher_results="${output_prefix}_philosopher.zip"
  }

  runtime 
  {
    docker : "broadcptac/pgdac_ms_philosopher:1"
    memory: "${if defined(ram_gb) then ram_gb else '26'}GB"
    disks : "local-disk ${if defined(local_disk_gb) then local_disk_gb else '120'} HDD"
    preemptible : "${if defined(num_preemptions) then num_preemptions else '0'}"
    cpu: "${if defined(num_cpus) then num_cpus else '4'}"
  }

  meta
  {
    author: "Felipe da Veiga Leprevost"
    email : "felipe@leprevost.com.br"
  }
}

task pgdac_tmt_integrator 
{
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

  output 
  {
    File results="results.zip"
  }

  runtime 
  {
    docker      : "broadcptac/pgdac_tmt_integrator:1"
    memory      : "${if defined(ram_gb) then ram_gb else '60'}GB"
    disks       : "local-disk ${if defined(local_disk_gb) then local_disk_gb else '100'} HDD"
    preemptible : "${if defined(num_preemptions) then num_preemptions else '0'}"
    cpu         : "${if defined(num_cpus) then num_cpus else '16'}"
  }

  meta 
  {
    author: "Felipe da Veiga Leprevost"
    email : "felipe@leprevost.com.br"
  }
}

workflow pgdac_philosopher_pipeline_workflow
{ 
  Boolean package
  String null_file="gs://broad-institute-gdac/GDAC_FC_NULL"
  File? database="gs://fc-secure-64796cb7-2cb9-42e2-a42a-741460a7efad/2019-10-31-td-RefSeq.20180629_Human_ucsc_hg38_cpdbnr_mito_264contams.fasta"
  File? fragger_params="gs://fc-secure-64796cb7-2cb9-42e2-a42a-741460a7efad/fragger.params"
  File? philosopher_params="gs://fc-secure-64796cb7-2cb9-42e2-a42a-741460a7efad/philosopher.yml"
  File? tmt_integrator_params="gs://fc-secure-64796cb7-2cb9-42e2-a42a-741460a7efad/tmt-i_param_v1.1.0.yml"

  File file_of_files
  Array[File] array_of_files = read_lines(file_of_files)
  
  call pgdac_msfragger
  {
    input: 
      package=package,
      package_name="msfragger-package",
      null_file="gs://broad-institute-gdac/GDAC_FC_NULL",
      mzML_files=array_of_files,
      params=fragger_params,
      database=database
  }

  #pepXML_archive=this.fragger_pepXML
  #philsosopher_pipeline_pkg=this.fragger_package

  call pgdac_philosopher 
  {
    input: 
      package=package,
      null_file=null_file,
      package_archive=pgdac_msfragger.philosopher_pipeline_pkg,
      mzML_files=array_of_files,
      pepXML_archive=pgdac_msfragger.pepXML_archive,
      database=database,
      params=philosopher_params
  }

  #philosopher_pipeline_pkg=this.philo_package
  #philosopher_results=this.philo_results

  call pgdac_tmt_integrator 
  {
    input: 
      database=database,
      philosopher_output=pgdac_philosopher.philosopher_results,
      params=tmt_integrator_params,
      database=database
  }
  
  output
  {
    File tmt_results=pgdac_tmt_integrator.results
  }
}
