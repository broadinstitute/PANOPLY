task panoply_philosopher_pipeline
{
  File philosopher_params
  File database
  Array[File] array_of_files
  
  String? output_prefix="philosopher"
  Boolean? package=true
  String? null_file
  Float? ram_gb = 128
  Int? local_disk_gb = 500
  Int? num_preemptions
  Int? num_cpus = 64

  command <<<
    . /etc/profile

    set -euo pipefail

    awk '{ gsub( "protein_database:", "protein_database: ${database}")}1' \
      ${philosopher_params} > philosopher.yml

    rootdir=$(dirname "${select_first( array_of_files )}")

    # Symlink project directory
    projdir=${output_prefix}_philosopher
    cp -Rs $rootdir $projdir

    # Go to project directory for generating annotation files and running
    # Philosopher
    mv philosopher.yml $projdir/philosopher.yml
    cd $projdir
    unzip *.zip

    # Run Philosopher
    philosopher pipeline --config philosopher.yml */

    # Run IonQuant
    java -Xmx32g -jar /IonQuant-1.7.17/IonQuant-1.7.17.jar --specdir $projdir/HeLa_subset $projdir/HeLa_subset/*.pepXML --psm $projdir/HeLa_subset/psm.tsv

    # Archive results
    cd -
    zip -r $projdir $projdir -x \*.d

    if ${package}; then
      package.sh -x broad-institute-gdac/\* \
                 -x \*.d \
                 -x $projdir.zip \
                 "philosopher_results.zip"
    fi

  >>>

  output 
  {
    File philosopher_pipeline_pkg="${if package then 'philosopher_results.zip' else null_file}"
    File philosopher_results="${output_prefix}_philosopher.zip"
  }

  runtime 
  {
    docker : "broadcptacdev/panoply_philosopher_pipeline:latest"
    memory: "${if defined(ram_gb) then ram_gb else '24'}GB"
    disks : "local-disk ${if defined(local_disk_gb) then local_disk_gb else '120'} HDD"
    preemptible : "${if defined(num_preemptions) then num_preemptions else '0'}"
    cpu: "${if defined(num_cpus) then num_cpus else '16'}"
  }

  meta
  {
    author: "Khoi Pham Munchic"
    email : "kpham@broadinstitute.org"
  }
}

workflow panoply_philosopher_pipeline_workflow
{ 
  File philosopher_params="gs://fc-secure-549d204e-cf25-4fef-898a-dc69ae22d1dd/philosopher_timsTOF.yml"
  File database="gs://fc-secure-549d204e-cf25-4fef-898a-dc69ae22d1dd/2019-10-31-td-RefSeq.20180629_Human_ucsc_hg38_cpdbnr_mito_264contams.fasta"
  File file_of_files="gs://fc-secure-549d204e-cf25-4fef-898a-dc69ae22d1dd/timsTOF_data/Proteome_d_files.txt"
  Array[File] array_of_files=read_lines(file_of_files)
  
  Boolean? package
  String? null_file="gs://broad-institute-gdac/GDAC_FC_NULL"

  call panoply_philosopher_pipeline 
  {
    input: 
      package=package,
      null_file=null_file,
      array_of_files=array_of_files,
      database=database,
      philosopher_params=philosopher_params
  }
}
