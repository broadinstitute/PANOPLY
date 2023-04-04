# DIA-NN spectral DIA search of raw files
**Version**: DIA-NN 1.8.1

## Description
This workflow wraps the DIA raw file search command in DIA-NN using library-free, combined library, or just empirical library (e.g., created by FragPipe). DIA-NN supports Thermo and Bruker files, and this workflow additionally automatically converts Thermo .raw into .mzML and timsTOF .d into .d.dia for running in Linux and for parallel processing.

A spectral library must be created prior to running this workflow (even if it's library-free library) -- either with `panoply_diann_speclib_gen` or directly in DIA-NN or FragPipe (via EasyPQP).

Step-by-step slides for running this workflow are available [here](https://docs.google.com/presentation/d/1nFFQ7DEMbGea4onjoXvKmemo2rVI4BDorlnt1igk31k/edit?usp=sharing).

## `panoply_diann_search` — workflow to search raw Thermo and Bruker files (.raw, .mzML, .d, .d.dia)
### Inputs
**Workflow inputs**
- `dia_files_folder` (Directory): Google Bucket path to folder containing raw DIA files to be searched (ex: "gs://fc-3f59fceb-8ce2-4855-9198-d6f6527cd8af/Experiment_1/")
- `speclib` (File): spectral library - library free, combined, or DDA library (.speclib format)
- `is_timsTOF` (Boolean): specify if data is from Bruker timsTOF (true) or Thermo Exploris (false)
- `additional_options` (String, optional): string specifying additional DIA-NN flags, refer to specific flag definitions [here](https://github.com/vdemichev/DiaNN#command-line-reference)

**DIA-NN parameters**
- `fragment_mz_max` (Integer, default=1800): maximum fragment m/z for DIA search
- `fragment_mz_min` (Integer, default=200): minimum fragment m/z for DIA search
- `missed_cleavages` (Integer, default=1): maximum number of missed cleavages
- `mod_C_carbamidomethylation` (Boolean, default=false): includes cysteine carbamidomethylation in DIA search 
- `mod_K_ubiquitination` (Boolean, default=false): includes lysine ubiquitination in DIA search 
- `mod_M_oxidation` (Boolean, default=false): includes methionine oxidation in DIA search 
- `mod_Nterm_acetylation` (Boolean, default=false): includes N-terminal acetylation in DIA search 
- `mod_Nterm_M_excision` (Boolean, default=false): removes methionine from N-terminal in DIA search 
- `mod_STY_phosphorylation` (Boolean, default=false): includes phosphorylation of serine, threonine, tyrosine in DIA search 
- `num_var_mods` (Integer, default=1): maximum number of variable modifications
- `peptide_len_max` (Integer, default=30): maximum precursor length for DIA search
- `peptide_len_min` (Integer, default=7): minimum precursor length for DIA search
- `precursor_charge_max` (Integer, default=4): maximum precursor charge for DIA search
- `precursor_charge_min` (Integer, default=1): minimum precursor charge for DIA search
- `precursor_mz_max` (Integer, default=1800): maximum precursor m/z for DIA search
- `precursor_mz_min` (Integer, default=300): minimum precursor m/z for DIA search
- To specify modifications not included by default, in `additional_options`, add flags of format `--var-mod UniMid:<NUMBER>,<MASS_SHIFT>,<RESIDUE>` (e.g., `--var-mod UniMod:259,8.014199,K  --var-mod UniMod:267,10.008269,R`)

**Skyline report (optional)**
- `skyline_generate_report` (Boolean, default=false): generate Skyline report from searched files
- `skyline_use_explicit_peak_bounds` (Boolean, default=true): if (`generate_skyline_report == true`), use the retention time apex and pick the peak boundaries (most users don't need to change this)
- `skyline_background_proteome_fasta` (File): Google Bucket path to background proteome file (.fasta)
- `skyline_template_zip` (File, optional): Google Bucket path to empty Skyline template (most users don't need to change this)

**Terra parameters**
- `local_disk_gb` (Integer, default=1000): GB of storage space in the compute instance
- `num_cpus` (Integer, default=32): number of cores in the compute instance
- `num_preemptions` (Integer, default=0): number 
- `ram_gb` (Integer, default=128): GB of RAM memory in the compute instance

### Outputs
From `dia_nn_first_pass`: contains search results for each individual run
- Each shard contains individual search of each run to enable parallel searching. The outputs are contained in `glob` folders;
    - `.quant` file for each run used for reanalysis with match-between-runs
    - Report tables for each run (aggregated in `diann_first_pass_out.zip` created by the final step, `dia_nn_match_between_runs`)

From `dia_nn_match_between_runs`:
- `diann_first_pass_out.zip`: contains search results and spectral library for individual runs
    - Long-format quality metrics table for each precursor in the individual run (`report.tsv`)
    - Quant table for each precursor ion in the individual run (`report.pr_matrix.tsv`)
    - Spectral library generated from that run (`.speclib`)
- `diann_out.zip`: contains final search results (after match-between-runs) and spectral library
    - Long-format quality metrics table for each precursor and run (`report.tsv`)
    - Quant table for each precursor ion across runs (`report.pr_matrix.tsv`)
    - Quant table for protein groups across runs (`report.pg_matrix.tsv`)
    - Quant table for gene groups across runs (`report.gg_matrix.tsv`)
    - Resulting spectral library of the entire search (`.speclib`)
- `dia_nn_match_between_runs.log`: Log file for DIA-NN search 

## Developer guide 
### Flow
1. `is_timsTOF` flag specifies how to convert files if needed
    - `true` would call conversion from timsTOF `.d` (a folder) to DIA-NN internal format `.d.dia` (a file)
    - `false` would call (parallelized) conversion from Exploris `.raw` to `.mzML` using [ThermoRawFileParser](https://github.com/compomics/ThermoRawFileParser)
2. Converted files are fed into first pass via scatter
    - Each first pass call receives a single `dia_file` and `speclib` to search with
    - Search parameters are passed down from main workflow
    - All generated `.quant` files and reports are gathered to pass to final step
3. First pass search results are fed into match-between-runs step
    - Converted and `.quant` files are used to do MBR (`--reanalyse` flag) 
    - Search parameters are passed down from main workflow
    - Final outputs and first pass outputs are stored

### Paths
- Working directory is `dia_nn`
- Input data will be in `dia_nn/data`
- Outputs will be in `dia_nn/out`

## References
- Demichev, V., Messner, C.B., Vernardis, S.I. et al. DIA-NN: neural networks and interference correction enable deep proteome coverage in high throughput. Nat Methods 17, 41–44 (2020). https://doi.org/10.1038/s41592-019-0638-x
- Steger, M., Demichev, V., Backman, M. et al. Time-resolved in vivo ubiquitinome profiling by DIA-MS reveals USP7 targets on a proteome-wide scale. Nat Commun 12, 5399 (2021). https://doi.org/10.1038/s41467-021-25454-1
- Demichev, V., Szyrwiel, L., Yu, F. et al. dia-PASEF data analysis using FragPipe and DIA-NN for deep proteomics of low sample amounts. Nat Commun 13, 3944 (2022). https://doi.org/10.1038/s41467-022-31492-0