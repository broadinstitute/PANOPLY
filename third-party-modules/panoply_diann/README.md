# DIA-NN spectral library generation and DIA search of raw files
**Version**: DIA-NN 1.8.1

## `panoply_diann_speclib_gen` — workflow to create a generated library or combined with DDA library
### Inputs
- `database` (File): Google Bucket path to background proteome file (.fasta)
- `speclib` (File, optional): Google Bucket path to secondary spectral library if you are making a combined library (.speclib)
- `additional_options` (String, optional): string specifying additional DIA-NN flags
### Outputs
- `speclib_predicted` (File): generated spectral library (library-free or combined)
- `speclib_prosit` (File): generated spectral library in PROSIT format
- `speclib_log` (File): logs during spectral library generation


## `panoply_diann_search` — workflow to search raw Thermo and Bruker files (.raw, .mzML, .d, .d.dia)
### Inputs
- `dia_files_folder` (Directory): Google Bucket path to folder containing raw DIA files to be searched (ex: "gs://fc-3f59fceb-8ce2-4855-9198-d6f6527cd8af/Experiment_1/")
- `speclib` (File): spectral library - generated, combined, or DDA (.speclib)
- `is_timsTOF` (Boolean): specify if data is from Bruker timsTOF (true) or Thermo Exploris (false)
- `additional_options` (String, optional): string specifying additional DIA-NN flags
### Outputs
From `dia_nn_match_between_runs`:
- diann_first_pass_out (File): zip containing search results and spectral library for individual runs
- diann_out (File): zip containing final search results (after match-between-runs) and spectral library

## Developer guide 
### Flow (`panoply_diann_speclib_gen`)
1. Only FASTA (and DDA library) are supplied with no data
2. Digest parameters and PTMs are specified
3. If no DDA library supplied, `--lib` flag will be ignored

### Flow (`panoply_diann_search`)
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
- Working directory is `/root/fragpipe`
- Input data will be in `/root/fragpipe/data`
- Outputs will be in `/root/fragpipe/out`
