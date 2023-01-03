version development
workflow panoply_dia_nn_search {
  input {
    File file_of_files  # must be .raw or .mzML (for Thermo) and .dia (for Bruker); .d folders cannot be specified
    File speclib
    Boolean is_timsTOF
    Array[File] files = read_lines(file_of_files)
    
    File? speclib_secondary
    File? config

    Int fragment_mz_min=300
    Int fragment_mz_max=1400
    Int precursor_mz_min=300
    Int precursor_mz_max=1400
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

    Boolean generate_spectral_library=false
    Boolean fasta_digest_library_free=false
    Boolean predict_spectra_rt_im=false
    String additional_options=""

    Int num_cpus=32
    Int ram_gb=128
    Int local_disk_gb=1000
    Int num_preemptions=0
  }
  
  if (!is_timsTOF) {
    scatter (file in files) {
      call convert_raw_to_mzml {
        input:
          file=file,

          num_cpus=4,
          ram_gb=16,
          local_disk_gb=local_disk_gb,
          num_preemptions=num_preemptions
      }
    }
  }
  
  scatter (file in (if is_timsTOF then files else select_first([convert_raw_to_mzml.converted_files]))) {
    call dia_nn_first_pass {
  	input:
      dia_file=file,
      speclib=speclib,
      speclib_secondary=speclib_secondary,
      config=config,

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

      generate_spectral_library=generate_spectral_library,
      fasta_digest_library_free=fasta_digest_library_free,
      predict_spectra_rt_im=predict_spectra_rt_im,
      additional_options=additional_options,

      num_cpus=num_cpus,
      ram_gb=ram_gb,
      local_disk_gb=local_disk_gb,
      num_preemptions=num_preemptions
    }
  }

  call dia_nn_match_between_runs {
    input:
      converted_files=(if is_timsTOF then files else select_first([convert_raw_to_mzml.converted_files])),
      quant_files=dia_nn_first_pass.quant,
      first_pass_outputs=dia_nn_first_pass.first_pass_out,

      speclib=speclib,
      speclib_secondary=speclib_secondary,
      config=config,

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

      generate_spectral_library=generate_spectral_library,
      fasta_digest_library_free=fasta_digest_library_free,
      predict_spectra_rt_im=predict_spectra_rt_im,
      additional_options=additional_options,

      num_cpus=num_cpus,
      ram_gb=ram_gb,
      local_disk_gb=local_disk_gb,
      num_preemptions=num_preemptions
  }
}

task dia_nn_first_pass {
  input {
    File dia_file
    File speclib
    File? speclib_secondary
    File? config
    
    Int fragment_mz_min=300
    Int fragment_mz_max=1400
    Int precursor_mz_min=300
    Int precursor_mz_max=1400
    Int precursor_charge_min=1
    Int precursor_charge_max=4
    Int missed_cleavages=1
    Int peptide_len_min=7
    Int peptide_len_max=30
    Int num_var_mods=1

    Boolean? mod_Nterm_M_excision
    Boolean? mod_C_carbamidomethylation
    Boolean? mod_M_oxidation
    Boolean? mod_STY_phosphorylation
    Boolean? mod_Nterm_acetylation
    Boolean? mod_K_ubiquitination

    Boolean? generate_spectral_library
    Boolean? fasta_digest_library_free
    Boolean? predict_spectra_rt_im
    String? additional_options

    Int? num_cpus
    Int? ram_gb
    Int? local_disk_gb
    Int? num_preemptions
  }
   
  command {
    . /etc/profile
    set -euo pipefail
    
    # copy from root to project dir
    projdir=dia_nn
    mkdir -p $projdir/data
    mkdir -p $projdir/out
    cp ${dia_file} $projdir/data/
    cp -s ${speclib} $projdir/

    dia=$(basename ${dia_file})
    lib=$(basename ${speclib})
    if [ -n "${speclib_secondary}" ]; then cp -s ${speclib_secondary} $projdir/; lib_second="--lib $(basename ${speclib_secondary})"; else lib_second=""; fi
    if [ -n "${config}" ]; then cp -s ${config} $projdir/; cfg="--cfg $(basename ${config})"; else cfg=""; fi

    # Parse flags
    if [ "${mod_Nterm_M_excision}" = true ]; then met_excision="--met-excision"; else met_excision=""; fi
    if [ "${mod_C_carbamidomethylation}" = true ]; then unimod4="--unimod4"; else unimod4=""; fi
    if [ "${mod_M_oxidation}" = true ]; then unimod35="--var-mod UniMod:35,15.994915,M"; else unimod35=""; fi
    if [ "${mod_Nterm_acetylation}" = true ]; then unimod1="--var-mod UniMod:1,42.010565,*n"; else unimod1=""; fi
    if [ "${mod_STY_phosphorylation}" = true ]; then unimod21="--var-mod UniMod:21,79.966331,STY"; else unimod21=""; fi
    if [ "${mod_K_ubiquitination}" = true ]; then unimod121="--var-mod UniMod:121,114.042927,K"; else unimod121=""; fi

    if [ "${generate_spectral_library}" = true ]; then gen_spec_lib="--gen-spec-lib"; else gen_spec_lib=""; fi
    if [ "${fasta_digest_library_free}" = true ]; then fasta_search="--fasta-search"; else fasta_search=""; fi
    if [ "${predict_spectra_rt_im}" = true ]; then predictor="--predictor"; else predictor=""; fi

    cd $projdir
    /usr/diann/1.8.1/diann-1.8.1 $cfg --dir "data" --lib $lib $lib_second --threads ${num_cpus} --verbose 1 \
      --out "out/report.tsv" --out-lib "out/spect_lib.tsv" --qvalue 0.01 --matrices --temp "data" $predictor --prosit --rt-profiling  \
      --min-fr-mz ${fragment_mz_min} --max-fr-mz ${fragment_mz_max} $met_excision --cut K*,R* --missed-cleavages ${missed_cleavages} \
      --min-pep-len ${peptide_len_min} --max-pep-len ${peptide_len_max} --min-pr-mz ${precursor_mz_min} --max-pr-mz ${precursor_mz_max} --min-pr-charge ${precursor_charge_min} --max-pr-charge ${precursor_charge_max} \
      $unimod4 --var-mods ${num_var_mods} $unimod35 $unimod1 $unimod21 $unimod121 \
      --smart-profiling --peak-center --no-ifs-removal --global-norm ${additional_options}

    # process output
    cd -
    mv -v $projdir/data/* .

    out_name=$(basename ${dia_file})-out
    cd $projdir && mv out $out_name && zip -r $out_name.zip $out_name && cd - && mv -v $projdir/$out_name.zip . 
  }

  output {
    File quant=select_first(glob("*.quant"))
    File first_pass_out=select_first(glob("*.zip"))
  }

  runtime {
    docker      : "broadcptacdev/panoply_dia_nn:latest"
    memory      : "${if defined(ram_gb) then ram_gb else '128'}GB"
    bootDiskSizeGb: 512
    disks       : "local-disk ${if defined(local_disk_gb) then local_disk_gb else '2000'} HDD"
    preemptible : "${if defined(num_preemptions) then num_preemptions else '0'}"
    cpu         : "${if defined(num_cpus) then num_cpus else '32'}"
  }

  meta {
    author: "Khoi Pham Munchic"
    email : "kpham@broadinstitute.org"
  }
}

task dia_nn_match_between_runs {
  input {
    Array[File] converted_files
    Array[File] quant_files
    Array[File] first_pass_outputs

    File speclib
    File? speclib_secondary
    File? config

    Int fragment_mz_min=300
    Int fragment_mz_max=1400
    Int precursor_mz_min=300
    Int precursor_mz_max=1400
    Int precursor_charge_min=1
    Int precursor_charge_max=4
    Int missed_cleavages=1
    Int peptide_len_min=7
    Int peptide_len_max=30
    Int num_var_mods=1

    Boolean? mod_Nterm_M_excision
    Boolean? mod_C_carbamidomethylation
    Boolean? mod_M_oxidation
    Boolean? mod_STY_phosphorylation
    Boolean? mod_Nterm_acetylation
    Boolean? mod_K_ubiquitination

    Boolean? generate_spectral_library
    Boolean? fasta_digest_library_free
    Boolean? predict_spectra_rt_im
    String? additional_options

    Int? num_cpus
    Int? ram_gb
    Int? local_disk_gb
    Int? num_preemptions
  }
   
  command {
    . /etc/profile
    set -euo pipefail
    
    # copy from root to project dir
    projdir=dia_nn
    mkdir -p $projdir/data
    mkdir -p $projdir/first_pass_out
    mkdir -p $projdir/out

    cp -s ~{sep(" ", converted_files)} $projdir/data/
    cp -s ~{sep(" ", quant_files)} $projdir/data/
    cp -s ${speclib} $projdir/
    lib=$(basename ${speclib})
    if [ -n "${speclib_secondary}" ]; then cp -s ${speclib_secondary} $projdir/; lib_second="--lib $(basename ${speclib_secondary})"; else lib_second=""; fi
    if [ -n "${config}" ]; then cp -s ${config} $projdir/; config="--cfg $(basename ${config})"; else config=""; fi

    cp -s  ~{sep(" ", first_pass_outputs)} $projdir/first_pass_out/
    cd $projdir/first_pass_out && unzip "*.zip" && rm *.zip

    # Parse flags
    if [ "${mod_Nterm_M_excision}" = true ]; then met_excision="--met-excision"; else met_excision=""; fi
    if [ "${mod_C_carbamidomethylation}" = true ]; then unimod4="--unimod4"; else unimod4=""; fi
    if [ "${mod_M_oxidation}" = true ]; then unimod35="--var-mod UniMod:35,15.994915,M"; else unimod35=""; fi
    if [ "${mod_Nterm_acetylation}" = true ]; then unimod1="--var-mod UniMod:1,42.010565,*n"; else unimod1=""; fi
    if [ "${mod_STY_phosphorylation}" = true ]; then unimod21="--var-mod UniMod:21,79.966331,STY"; else unimod21=""; fi
    if [ "${mod_K_ubiquitination}" = true ]; then unimod121="--var-mod UniMod:121,114.042927,K"; else unimod121=""; fi


    if [ "${generate_spectral_library}" = true ]; then gen_spec_lib="--gen-spec-lib"; else gen_spec_lib=""; fi
    if [ "${fasta_digest_library_free}" = true ]; then fasta_search="--fasta-search"; else fasta_search=""; fi
    if [ "${predict_spectra_rt_im}" = true ]; then predictor="--predictor"; else predictor=""; fi

    # Run DIA-NN
    cd .. # we're in projdir
    /usr/diann/1.8.1/diann-1.8.1 $config --dir "data" --lib $lib $lib_second --threads ${num_cpus} --verbose 1 \
      --out "out/report.tsv" --out-lib "out/spect_lib.tsv" --qvalue 0.01 --matrices --temp "data" $predictor --prosit --rt-profiling \
      --min-fr-mz ${fragment_mz_min} --max-fr-mz ${fragment_mz_max} $met_excision --cut K*,R* --missed-cleavages ${missed_cleavages} \
      --min-pep-len ${peptide_len_min} --max-pep-len ${peptide_len_max} --min-pr-mz ${precursor_mz_min} --max-pr-mz ${precursor_mz_max} --min-pr-charge ${precursor_charge_min} --max-pr-charge ${precursor_charge_max} \
      $unimod4 --var-mods ${num_var_mods} $unimod35 $unimod1 $unimod21 $unimod121 \
      --use-quant --reanalyse --smart-profiling --peak-center --no-ifs-removal --global-norm ${additional_options}

    # Archive results
    cd ..
    zip -r "diann_first_pass_out.zip" $projdir/first_pass_out
    zip -r "diann_out.zip" $projdir/out
  }

  output {
    File diann_first_pass_out="diann_first_pass_out.zip"
    File diann_out="diann_out.zip"
  }

  runtime {
    docker      : "broadcptacdev/panoply_dia_nn:latest"
    memory      : "${if defined(ram_gb) then ram_gb else '128'}GB"
    bootDiskSizeGb: 512
    disks       : "local-disk ${if defined(local_disk_gb) then local_disk_gb else '2000'} HDD"
    preemptible : "${if defined(num_preemptions) then num_preemptions else '0'}"
    cpu         : "${if defined(num_cpus) then num_cpus else '32'}"
  }

  meta {
    author: "Khoi Pham Munchic"
    email : "kpham@broadinstitute.org"
  }
}

task convert_raw_to_mzml {
  input {
    File file

    Int? num_cpus
    Int? ram_gb
    Int? local_disk_gb
    Int? num_preemptions
  }
   
  command {
    . /etc/profile
    set -euo pipefail
    
    # copy from root to project dir
    projdir=dia_nn
    mkdir -p $projdir/data

    cd $projdir
    cp ~{file} data/

    if ls data/*.raw &> /dev/null; then
      mono /ThermoRawFileParser/ThermoRawFileParser.exe -d data/ -o data/
    else
      echo "No .raw file is present. Skipping .raw to .mzML conversion."
    fi

    cd ..
    find $projdir/data -name '*.mzML' | xargs mv -t .
  }

  output {
    File converted_files = select_first(glob("*.mzML"))
  }

  runtime {
    docker      : "broadcptacdev/panoply_dia_nn:latest"
    memory      : "${if defined(ram_gb) then ram_gb else '16'}GB"
    bootDiskSizeGb: 128
    disks       : "local-disk ${if defined(local_disk_gb) then local_disk_gb else '750'} HDD"
    preemptible : "${if defined(num_preemptions) then num_preemptions else '0'}"
    cpu         : "${if defined(num_cpus) then num_cpus else '4'}"
  }

  meta {
    author: "Khoi Pham Munchic"
    email : "kpham@broadinstitute.org"
  }
}
