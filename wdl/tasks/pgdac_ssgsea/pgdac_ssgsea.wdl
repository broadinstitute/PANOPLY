task pgdac_ssgsea {

	File input_ds
	
  ## parameters to create gene-centric or single-site-centric 
  ## GCT files for ssGSEA / PTM-SEA
	String level="ssc"
  String id_type="sm"
  String id_type_out="uniprot"
  String acc_type="refseq"
  String seqwin_col="VMsiteFlanks"
  String gene_col="geneSymbol"
  String SGT_col="subgroupNum"
  Boolean loc=true
  String mode="median"
  String mod_res="S|T|Y"
  String mod_type="p"
  Boolean use_as_is=false
	
	## ssGSEA / PTM-SEA parameters below
	File gene_set_database
	String output_prefix
	String sample_norm_type="rank"
	String correl_type="z.score"
	String statistic="area.under.RES"
	String output_score_type="NES"
	Float weight=0.75
	Int min_overlap=5
	Int nperm=1000
	Boolean global_fdr=true

  ## VM parameters
	Int memory
	Int disk_space
	Int num_threads
	
	command {
		set -euo pipefail
		
		# prepare GCT file
		/home/pgdac/src/preprocessGCT.R -i ${input_ds} -l ${level} -t ${id_type} -o ${id_type_out} -a ${acc_type} -s ${seqwin_col} --gene_symbol ${gene_col} -v ${SGT_col} -d ${loc} -m ${mode} -r "${mod_res}" -p '${mod_type}' -u ${use_as_is} -z /home/pgdac/src
		
		# update path to input_ds
		input_ds=`cat fn.out`
		
		# run ssgsea/ptm-sea
		/home/pgdac/ssgsea-cli.R -i ${input_ds} -d ${gene_set_database} -o ${output_prefix} -n ${sample_norm_type} -w ${weight} -c ${correl_type} -t ${statistic} -s ${output_score_type} -p ${nperm} -m ${min_overlap} -g ${global_fdr}
		
		# archive result files
		find * -regextype posix-extended -regex "^signature_gct/.*.gct$|^${output_prefix}.*.gct$|^.*.log.txt$|^.*parameters.txt$" -print0 | tar -czvf ${output_prefix}.tar.gz --null -T -
		}

	output {
		# Outputs defined here
		File results="${output_prefix}.tar.gz"
		}

	runtime {
		docker : "broadcptac/pgdac_ssgsea:5"
		memory : "${memory}GB"
		disks : "local-disk ${disk_space} HDD"
		cpu : "${num_threads}"
	}

	meta {
		author : "Karsten Krug"
		email : "karsten@broadinstitute.org"
	}

}

workflow pgdac_ssgsea_workflow {
	call pgdac_ssgsea
}
