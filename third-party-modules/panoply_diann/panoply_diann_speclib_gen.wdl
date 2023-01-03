task diann_speclib_gen
{
  File database
  File? speclib

  Int? fragment_mz_min=300
  Int? fragment_mz_max=1400
  Int? precursor_mz_min=300
  Int? precursor_mz_max=1400
  Int? precursor_charge_min=1
  Int? precursor_charge_max=4
  Int? missed_cleavages=1
  Int? peptide_len_min=7
  Int? peptide_len_max=30
  Int? num_var_mods=1

  Boolean? mod_Nterm_M_excision=true
  Boolean? mod_C_carbamidomethylation=true
  Boolean? mod_M_oxidation=false
  Boolean? mod_STY_phosphorylation=false
  Boolean? mod_Nterm_acetylation=false

  Boolean? generate_spectral_library=true
  Boolean? fasta_digest_library_free=true
  Boolean? predict_spectra_rt_im=true
  String? additional_options=""

  Int? num_cpus
  Int? ram_gb
  Int? local_disk_gb
  Int? num_preemptions
   
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
    if [ "${mod_Nterm_acetylation}" = true ]; then unimod1="--var-mod UniMod:1,42.010565,*n"; else unimod1=""; fi
    if [ "${mod_STY_phosphorylation}" = true ]; then unimod21="--var-mod UniMod:21,79.966331,STY"; else unimod21=""; fi

    if [ "${generate_spectral_library}" = true ]; then gen_spec_lib="--gen-spec-lib"; else gen_spec_lib=""; fi
    if [ "${fasta_digest_library_free}" = true ]; then fasta_search="--fasta-search"; else fasta_search=""; fi
    if [ "${predict_spectra_rt_im}" = true ]; then predictor="--predictor"; else predictor=""; fi

    # Run DIA-NN
    cd $projdir
    /usr/diann/1.8.1/diann-1.8.1 --lib $lib --threads ${num_cpus} --verbose 1 \
      --out "out/report.tsv" --qvalue 0.01 --matrices --temp "out" --out-lib "out/speclib.tsv" $gen_spec_lib --prosit $predictor --fasta $db $fasta_search \
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
    memory      : "${if defined(ram_gb) then ram_gb else '128'}GB"
    disks       : "local-disk ${if defined(local_disk_gb) then local_disk_gb else '500'} HDD"
    preemptible : "${if defined(num_preemptions) then num_preemptions else '0'}"
    cpu         : "${if defined(num_cpus) then num_cpus else '32'}"
  }

  meta {
    author: "Khoi Pham Munchic"
    email : "kpham@broadinstitute.org"
  }
}

workflow panoply_diann_speclib_gen {
  Int? num_cpus=32
  Int? ram_gb=128
  Int? local_disk_gb=500
  Int? num_preemptions=0

  call diann_speclib_gen
  {
  	input:
      num_cpus=num_cpus,
      ram_gb=ram_gb,
      local_disk_gb=local_disk_gb,
      num_preemptions=num_preemptions
  }
}
