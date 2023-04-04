# DIA-NN spectral library generation
**Version**: DIA-NN 1.8.1

## Description
This workflow wraps the spectral library prediction command for DIA search in DIA-NN (deep learning-based spectra, RTs and IMs prediction). There are two modes of generation: library-free (in which a .FASTA proteome file is supplied) and combined (in which both proteome .FASTA and a DDA library is provided). For the combined library, DIA-NN supports `.csv`, `.tsv`, `.xls`, `.txt`, `.speclib`, `.sptxt`, `.msp` library files. The resulting library is generated in `.speclib` and PROSIT `.txt` formats.

Step-by-step slides for running this workflow are available [here](https://docs.google.com/presentation/d/1nFFQ7DEMbGea4onjoXvKmemo2rVI4BDorlnt1igk31k/edit?usp=sharing).

## `panoply_diann_speclib_gen` — workflow to create a generated library or combined with DDA library
### Inputs
**Workflow inputs**
- `database` (File): Google Bucket path to background proteome file (.fasta)
- `speclib` (File, optional): Google Bucket path to secondary spectral library if you are making a combined library (.speclib)
- `additional_options` (String, optional): string specifying additional DIA-NN flags, refer to specific flag definitions [here](https://github.com/vdemichev/DiaNN#command-line-reference)

**DIA-NN parameters**
- `fragment_mz_max` (Integer, default=1800): maximum fragment m/z for the in silico library generation
- `fragment_mz_min` (Integer, default=200): minimum fragment m/z for the in silico library generation
- `missed_cleavages` (Integer, default=1): maximum number of missed cleavages
- `mod_C_carbamidomethylation` (Boolean, default=false): includes cysteine carbamidomethylation in in silico library generation 
- `mod_K_ubiquitination` (Boolean, default=false): includes lysine ubiquitination in in silico library generation 
- `mod_M_oxidation` (Boolean, default=false): includes methionine oxidation in in silico library generation 
- `mod_Nterm_acetylation` (Boolean, default=false): includes N-terminal acetylation in in silico library generation 
- `mod_Nterm_M_excision` (Boolean, default=false): removes methionine from N-terminal in in silico library generation 
- `mod_STY_phosphorylation` (Boolean, default=false): includes phosphorylation of serine, threonine, tyrosine in in silico library generation 
- `num_var_mods` (Integer, default=1): maximum number of variable modifications
- `peptide_len_max` (Integer, default=30): maximum precursor length for in silico library generation
- `peptide_len_min` (Integer, default=7): minimum precursor length for in silico library generation
- `precursor_charge_max` (Integer, default=4): maximum precursor charge for in silico library generation
- `precursor_charge_min` (Integer, default=1): minimum precursor charge for in silico library generation
- `precursor_mz_max` (Integer, default=1800): maximum precursor m/z for in silico library generation
- `precursor_mz_min` (Integer, default=300): minimum precursor m/z for in silico library generation
- To specify modifications not included by default, in `additional_options`, add flags of format `--var-mod UniMid:<NUMBER>,<MASS_SHIFT>,<RESIDUE>` (e.g., `--var-mod UniMod:259,8.014199,K  --var-mod UniMod:267,10.008269,R`)

**Terra parameters**
- `local_disk_gb` (Integer, default=1000): GB of storage space in the compute instance
- `num_cpus` (Integer, default=32): number of cores in the compute instance
- `num_preemptions` (Integer, default=0): number 
- `ram_gb` (Integer, default=128): GB of RAM memory in the compute instance

### Outputs
- `speclib.predicted.speclib` (File): generated spectral library (library-free or combined)
- `speclib.prosit.csv` (File): generated spectral library in PROSIT format
- `speclib.log.txt` (File): DIA-NN log during spectral library generation
- `dia_nn_speclib_gen.log`: Log file for DIA-NN library generation 

## Developer guide 
### Flow
1. Only FASTA (and DDA library, if combined library) are supplied with no data
2. Digest parameters and PTMs are specified
3. If no DDA library supplied, `--lib` flag will be ignored
4. DIA-NN command call without supplying files to search

### Paths
- Working directory is `dia_nn`
- Input data will be in `dia_nn/data`
- Outputs will be in `dia_nn/out`

## References
- Demichev, V., Messner, C.B., Vernardis, S.I. et al. DIA-NN: neural networks and interference correction enable deep proteome coverage in high throughput. Nat Methods 17, 41–44 (2020). https://doi.org/10.1038/s41592-019-0638-x
- Steger, M., Demichev, V., Backman, M. et al. Time-resolved in vivo ubiquitinome profiling by DIA-MS reveals USP7 targets on a proteome-wide scale. Nat Commun 12, 5399 (2021). https://doi.org/10.1038/s41467-021-25454-1
- Demichev, V., Szyrwiel, L., Yu, F. et al. dia-PASEF data analysis using FragPipe and DIA-NN for deep proteomics of low sample amounts. Nat Commun 13, 3944 (2022). https://doi.org/10.1038/s41467-022-31492-0