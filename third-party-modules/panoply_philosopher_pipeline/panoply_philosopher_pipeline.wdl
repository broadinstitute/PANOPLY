task panoply_philosopher_pipeline
{
  File philosopher_params
  File database
  File ms_data_zip
  
  String? output_prefix
  Boolean? package
  String? null_file
  Int? num_preemptions
  Int? num_cpus
  Int? ram_gb
  Int? local_disk_gb
  
  # MSFragger
  Int? MSFragger__num_threads = 24
  Int? MSFragger__memory = 64
  String? MSFragger__decoy_tag
  String? MSFragger__contaminant_tag
  Int? MSFragger__precursor_mass_lower
  Int? MSFragger__precursor_mass_upper
  String? MSFragger__variable_mod_01
  String? MSFragger__variable_mod_02
  String? MSFragger__variable_mod_03
  String? MSFragger__variable_mod_04
  String? MSFragger__variable_mod_05
  String? MSFragger__variable_mod_06
  String? MSFragger__variable_mod_07
  Int? MSFragger__max_variable_mods_per_peptide
  Int? MSFragger__digest_min_length
  Int? MSFragger__digest_max_length
  String? MSFragger__digest_mass_range

  # Filter
  Float? Filter__psmFDR
  Float? Filter__peptideFDR
  Float? Filter__ionFDR
  Float? Filter__proteinFDR

  command <<<
    . /etc/profile
    set -euo pipefail

    # Symlink project directory
    projdir=${output_prefix}_philosopher
    mkdir $projdir
    cp -s ${ms_data_zip} $projdir/

    # # Create updated philosopher parameters YAML
    # Checks if a variable is set then uses yq to edit the respective parameter in the yaml file
    # String and numeric are unwrapped differently: ex. "decoy_tag" as String and "ram_gb" as Int
    database=${database} yq e '.["Database Search"].protein_database = env(database)' ${philosopher_params} > philosopher.yml
    if [ -n "${MSFragger__memory}" ]; then memory=${MSFragger__memory} yq e '.["Database Search"].msfragger.memory = env(memory)' -i philosopher.yml; fi
    if [ -n "${MSFragger__num_threads}" ]; then num_threads=${MSFragger__num_threads} yq e '.["Database Search"].msfragger.num_threads = env(num_threads)' -i philosopher.yml; fi
    if [ -n "${MSFragger__decoy_tag}" ]; then decoy_tag="${MSFragger__decoy_tag}" yq e '.["Database Search"].decoy_tag = env(decoy_tag)' -i philosopher.yml; fi
    if [ -n "${MSFragger__contaminant_tag}" ]; then contaminant_tag="${MSFragger__contaminant_tag}" yq e '.["Database Search"].contaminant_tag = env(contaminant_tag)' -i philosopher.yml; fi

    if [ -n "${MSFragger__precursor_mass_lower}" ]; then precursor_mass_lower=${MSFragger__precursor_mass_lower} yq e '.["Database Search"].msfragger.precursor_mass_lower = env(precursor_mass_lower)' -i philosopher.yml; fi
    if [ -n "${MSFragger__precursor_mass_upper}" ]; then precursor_mass_upper=${MSFragger__precursor_mass_upper} yq e '.["Database Search"].msfragger.precursor_mass_upper = env(precursor_mass_upper)' -i philosopher.yml; fi
    if [ -n "${MSFragger__variable_mod_01}" ]; then variable_mod_01="${MSFragger__variable_mod_01}" yq e '.["Database Search"].msfragger.variable_mod_01 = env(variable_mod_01)' -i philosopher.yml; fi
    if [ -n "${MSFragger__variable_mod_02}" ]; then variable_mod_02="${MSFragger__variable_mod_02}" yq e '.["Database Search"].msfragger.variable_mod_02 = env(variable_mod_02)' -i philosopher.yml; fi
    if [ -n "${MSFragger__variable_mod_03}" ]; then variable_mod_03="${MSFragger__variable_mod_03}" yq e '.["Database Search"].msfragger.variable_mod_03 = env(variable_mod_03)' -i philosopher.yml; fi
    if [ -n "${MSFragger__variable_mod_04}" ]; then variable_mod_04="${MSFragger__variable_mod_04}" yq e '.["Database Search"].msfragger.variable_mod_04 = env(variable_mod_04)' -i philosopher.yml; fi
    if [ -n "${MSFragger__variable_mod_05}" ]; then variable_mod_05="${MSFragger__variable_mod_05}" yq e '.["Database Search"].msfragger.variable_mod_05 = env(variable_mod_05)' -i philosopher.yml; fi
    if [ -n "${MSFragger__variable_mod_06}" ]; then variable_mod_06="${MSFragger__variable_mod_06}" yq e '.["Database Search"].msfragger.variable_mod_06 = env(variable_mod_06)' -i philosopher.yml; fi
    if [ -n "${MSFragger__variable_mod_07}" ]; then variable_mod_07="${MSFragger__variable_mod_07}" yq e '.["Database Search"].msfragger.variable_mod_07 = env(variable_mod_07)' -i philosopher.yml; fi
    if [ -n "${MSFragger__max_variable_mods_per_peptide}" ]; then max_variable_mods_per_peptide=${MSFragger__max_variable_mods_per_peptide} yq e '.["Database Search"].msfragger.max_variable_mods_per_peptide = env(max_variable_mods_per_peptide)' -i philosopher.yml; fi
    if [ -n "${MSFragger__digest_min_length}" ]; then digest_min_length=${MSFragger__digest_min_length} yq e '.["Database Search"].msfragger.digest_min_length = env(digest_min_length)' -i philosopher.yml; fi
    if [ -n "${MSFragger__digest_max_length}" ]; then digest_max_length=${MSFragger__digest_max_length} yq e '.["Database Search"].msfragger.digest_max_length = env(digest_max_length)' -i philosopher.yml; fi
    if [ -n "${MSFragger__digest_mass_range}" ]; then digest_mass_range="${MSFragger__digest_mass_range}" yq e '.["Database Search"].msfragger.digest_mass_range = env(digest_mass_range)' -i philosopher.yml; fi

    if [ -n "${Filter__psmFDR}" ]; then psmFDR=${Filter__psmFDR} yq e '.["FDR Filtering"].psmFDR = env(psmFDR)' -i philosopher.yml; fi
    if [ -n "${Filter__peptideFDR}" ]; then peptideFDR=${Filter__peptideFDR} yq e '.["FDR Filtering"].peptideFDR = env(peptideFDR)' -i philosopher.yml; fi
    if [ -n "${Filter__ionFDR}" ]; then ionFDR=${Filter__ionFDR} yq e '.["FDR Filtering"].ionFDR = env(ionFDR)' -i philosopher.yml; fi
    if [ -n "${Filter__proteinFDR}" ]; then proteinFDR=${Filter__proteinFDR} yq e '.["FDR Filtering"].proteinFDR = env(proteinFDR)' -i philosopher.yml; fi

    # Go to project directory for generating annotation files and running Philosopher
    mv philosopher.yml $projdir/philosopher.yml
    cd $projdir
    unzip *.zip

    # Run Philosopher
    philosopher pipeline --config philosopher.yml */

    # Run IonQuant
    # java -Xmx32g -jar /IonQuant-1.7.17/IonQuant-1.7.17.jar --specdir $projdir/HeLa_subset $projdir/HeLa_subset/*.pepXML --psm $projdir/HeLa_subset/psm.tsv

    # Archive results
    cd -
    zip -r $projdir $projdir -x \*.zip
    zip -r "philosopher_results.zip" . 


  >>>

  output 
  {
    # File philosopher_yaml="philosopher.yml"
    File philosopher_results="${output_prefix}_philosopher.zip"
    File philosopher_pipeline_pkg="philosopher_results.zip"
  }

  runtime 
  {
    docker : "broadcptacdev/panoply_philosopher_pipeline:latest"
    memory: "${if defined(ram_gb) then ram_gb else '128'}GB"
    disks : "local-disk ${if defined(local_disk_gb) then local_disk_gb else '500'} HDD"
    preemptible : "${if defined(num_preemptions) then num_preemptions else '0'}"
    cpu: "${if defined(num_cpus) then num_cpus else '32'}"
  }

  meta
  {
    author: "Khoi Pham Munchic"
    email : "kpham@broadinstitute.org"
  }
}

workflow panoply_philosopher_pipeline_workflow
{ 
  File philosopher_params
  File database
  # File file_of_files
  # Array[File] array_of_files=read_lines(file_of_files)
  File ms_data_zip
  
  Boolean? package=true
  String? output_prefix="philosopher"
  String? null_file="gs://broad-institute-gdac/GDAC_FC_NULL"
  Int? num_preemptions
  Int? num_cpus=32
  Int? ram_gb=128
  Int? local_disk_gb=500

  call panoply_philosopher_pipeline 
  {
    input: 
      philosopher_params=philosopher_params,
      database=database,
      ms_data_zip=ms_data_zip,

      package=package,
      output_prefix=output_prefix,
      null_file=null_file,
      num_preemptions=num_preemptions,
      num_cpus=num_cpus,
      ram_gb=ram_gb,
      local_disk_gb=local_disk_gb
  }
}
