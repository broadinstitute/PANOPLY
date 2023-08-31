workflow panoply_diann_speclib_gen {
  File database
  File? speclib

  Int fragment_mz_min=200
  Int fragment_mz_max=1800
  Int precursor_mz_min=300
  Int precursor_mz_max=1800
  Int precursor_charge_min=1
  Int precursor_charge_max=4
  Int missed_cleavages=1
  Int peptide_len_min=7
  Int peptide_len_max=30
  Int num_var_mods=1

  Boolean mod_Nterm_M_excision=false
  Boolean mod_C_carbamidomethylation=false
  Boolean mod_M_oxidation=false
  Boolean mod_STY_phosphorylation=false
  Boolean mod_Nterm_acetylation=false
  Boolean mod_K_ubiquitination=false

  String additional_options=""

  Int num_cpus=32
  Int ram_gb=128
  Int local_disk_gb=500
  Int num_preemptions=0

  call diann_speclib_gen
  {
  	input:
      database=database,
      speclib=speclib,

      fragment_mz_min=fragment_mz_min,
      fragment_mz_max=fragment_mz_max,
      precursor_mz_min=precursor_mz_min,
      precursor_mz_max=precursor_mz_max,
      precursor_charge_min=precursor_charge_min,
      precursor_charge_max=precursor_charge_max,
      missed_cleavages=missed_cleavages,
      peptide_len_min=peptide_len_min,
      peptide_len_max=peptide_len_max,
      num_var_mods=num_var_mods,

      mod_Nterm_M_excision=mod_Nterm_M_excision,
      mod_C_carbamidomethylation=mod_C_carbamidomethylation,
      mod_M_oxidation=mod_M_oxidation,
      mod_STY_phosphorylation=mod_STY_phosphorylation,
      mod_Nterm_acetylation=mod_Nterm_acetylation,
      mod_K_ubiquitination=mod_K_ubiquitination,

      additional_options=additional_options,

      num_cpus=num_cpus,
      ram_gb=ram_gb,
      local_disk_gb=local_disk_gb,
      num_preemptions=num_preemptions
  }
}

task diann_speclib_gen
{
  File database
  File? speclib

  Int fragment_mz_min=200
  Int fragment_mz_max=1800
  Int precursor_mz_min=300
  Int precursor_mz_max=1800
  Int precursor_charge_min=1
  Int precursor_charge_max=4
  Int missed_cleavages=1
  Int peptide_len_min=7
  Int peptide_len_max=30
  Int num_var_mods=1

  Boolean mod_Nterm_M_excision
  Boolean mod_C_carbamidomethylation
  Boolean mod_M_oxidation
  Boolean mod_STY_phosphorylation
  Boolean mod_Nterm_acetylation
  Boolean mod_K_ubiquitination

  String? additional_options

  Int num_cpus=32
  Int ram_gb=128
  Int local_disk_gb=500
  Int num_preemptions=0
   
  command <<<
    . /etc/profile
    set -euo pipefail
    
    # copy from root to project dir
    projdir=dia_nn
    mkdir -p $projdir/data
    mkdir -p $projdir/out
    cp -s ${database} $projdir/

    db=$(basename ${database})
    if [ -n "${speclib}" ]
    then
      cp -s ${speclib} $projdir/
      lib=$(basename ${speclib})
    else
      lib=""
    fi

    # Parse flags
    if [ "${mod_Nterm_M_excision}" = true ]; then met_excision="--met-excision"; else met_excision=""; fi
    if [ "${mod_C_carbamidomethylation}" = true ]; then unimod4="--unimod4"; else unimod4=""; fi
    if [ "${mod_M_oxidation}" = true ]; then unimod35="--var-mod UniMod:35,15.994915,M"; else unimod35=""; fi
    if [ "${mod_STY_phosphorylation}" = true ]; then unimod21="--var-mod UniMod:21,79.966331,STY"; else unimod21=""; fi
    if [ "${mod_Nterm_acetylation}" = true ]; then unimod1="--var-mod UniMod:1,42.010565,*n"; else unimod1=""; fi
    if [ "${mod_K_ubiquitination}" = true ]; then unimod121="--var-mod UniMod:121,114.042927,K"; else unimod121=""; fi

    # Run DIA-NN
    cd $projdir
    /usr/diann/1.8.1/diann-1.8.1 --lib $lib --threads ${num_cpus} --verbose 1 \
      --out "out/report.tsv" --qvalue 0.01 --matrices --temp "out" --out-lib "out/speclib.tsv" --gen-spec-lib --prosit --predictor --fasta $db --fasta-search \
      --min-fr-mz ${fragment_mz_min} --max-fr-mz ${fragment_mz_max} $met_excision --cut K*,R* --missed-cleavages ${missed_cleavages} \
      --min-pep-len ${peptide_len_min} --max-pep-len ${peptide_len_max} --min-pr-mz ${precursor_mz_min} --max-pr-mz ${precursor_mz_max} --min-pr-charge ${precursor_charge_min} --max-pr-charge ${precursor_charge_max} \
      $unimod4 --var-mods ${num_var_mods} $unimod35 $unimod1 $unimod21 --monitor-mod UniMod:1 \
      --use-quant --reanalyse --smart-profiling --peak-center --no-ifs-removal --global-norm \
      ${additional_options}

    cd -
    mv $projdir/out/speclib.predicted.speclib .
    mv $projdir/out/speclib.prosit.csv .
    mv $projdir/out/speclib.log.txt .
  >>>

  output {
    File speclib_predicted="speclib.predicted.speclib"
    File speclib_prosit="speclib.prosit.csv"
    File speclib_log="speclib.log.txt"
  }

  runtime {
    docker      : "broadcptacdev/panoply_dia_nn:latest"
    memory      : "${ram_gb}GB"
    disks       : "local-disk ${local_disk_gb} HDD"
    preemptible : num_preemptions
    cpu         : num_cpus
  }

  meta {
    author: "Khoi Pham Munchic"
    email : "kpham@broadinstitute.org"
  }
}
