# ```panoply_quilts```

## Description

Workflow to perform QUILTS (https://github.com/ekawaler/pyQUILTS) on sample data (somatic, germline, splice junction, and fusion). QUILTS produces a custom protein database which can be used in mass spectrometry searches.

## Input

### Required inputs:

* ```junction_file_type```: (String default = 'star') String to indicate what junction finding software was used (only accepts star, MapSplice, or tophat).
* ```out_file```: (String default = 'quilts_results.tar.gz') Name of the output file produced.
* ```output_dir```: (String default = '/src/QUILTS/output') Where the QUILTS data will be deposited on the docker.
* ```reference_genome```: (File default = 'gs://fc-5a1bdedd-17aa-4a16-b661-bc16902fc4e0/genome/emily_ensembl_hg38.tar.gz') A compressed file containing the genome to be used by QUILTS. See the pyQUILTS repository for more info on making a genome for QUILTS use. By default an ensembl hg38 genome is used. More premade genomes are available on the Terra workspace 'QUILTS_data' google bucket. 
* ```reference_proteome```: (File default = 'gs://fc-5a1bdedd-17aa-4a16-b661-bc16902fc4e0/proteome/emily_ensembl_v100_hg38.tar.gz') A compressed file containing the proteome to be used by QUILTS. See the pyQUILTS repository for more info on making a proteome for QUILTS use. By default an ensembl version 100 hg38 proteome is used. More premade proteomes are available on the Terra workspace 'QUILTS_data' google bucket.
* At least one of:
* ```input_somatic_vcfs```: (array of `.vcf` files) Somatic files must be in `.vcf` format. This input can take an array of files but only input arrays of files if you wish for the output to be a merged file containing all results. If you wish to run QUILTS on a per sample basis please enter only one data file in these input array fields.
* ```input_germline_vcfs```: (array of `.vcf` files) Germline files must be in `.vcf` format. This input can take an array of files but only input arrays of files if you wish for the output to be a merged file containing all results. If you wish to run QUILTS on a per sample basis please enter only one data file in these input array fields.
* ```input_gene_fusions_files```: (array of `.txt` files) Gene fusions files must be a `.txt` in the following format:

FusionName	LeftBreakpoint	RightBreakpoint	Sample	JunctionReadCount
GENE1—GENE2	chr19:1000000:-	chr16:500003:+	Sample1	5 

This input can take an array of files but only input arrays of files if you wish for the output to be a merged file containing all results. If you wish to run QUILTS on a per sample basis please enter only one data file in these input array fields.
* ```input_splice_junctions_files```: (array of either `.txt`, `.tab`, or `.bed` files) Splice junction files must be either `.txt` for MapSplice results, `.tab` for STAR results, or `.bed` for tophat results. This input can take an array of files but only input arrays of files if you wish for the output to be a merged file containing all results. If you wish to run QUILTS on a per sample basis please enter only one data file in these input array fields.


### Optional inputs:

* ```threshB```: (integer) Please see the pyQUILTS repository for more information on this input.
* ```threshD```: (integer) Please see the pyQUILTS repository for more information on this input.
* ```threshN```: (integer) Please see the pyQUILTS repository for more information on this input.
* ```variant_quality_threshold```: (integer) Please see the pyQUILTS repository for more information on this input.

## Output

* ```outputs```: (File `.tar.gz`) The name of this output depends on the input to ```out_file``` above. This file contains the results of the QUILTS run in a compressed format.

