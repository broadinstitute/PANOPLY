task quilts {
	Array[File]? input_somatic_vcfs
	Array[File]? input_germline_vcfs
	Array[File]? input_splice_junctions_files # .bed .txt or .tab (see below)
	String? junction_file_type   # options: 'mapsplice', 'tophat', 'star' 
	Array[File]? input_gene_fusions_files     # .txt file of doc specified formatting
	String output_dir='/src/QUILTS/output' # default results dir set in the docker
	File reference_genome='gs://fc-5a1bdedd-17aa-4a16-b661-bc16902fc4e0/genome/emily_ensembl_hg38.tar.gz'   # set to hg38 by default
	File reference_proteome='gs://fc-5a1bdedd-17aa-4a16-b661-bc16902fc4e0/proteome/emily_ensembl_v100_hg38.tar.gz' # set to hg38 v100 (gencodev34 equivalent) by default
	String out_file='quilts_results.tar.gz' # default output name so the wdl can find the results

	Int? threshB                   # integer - minimum number of reads to support a splice junction with conserved exon boundaries
	Int? threshD				   # integer - minimum number of reads to support a splice junction with only the donor exon boundary conserved
	Int? threshN				   # integer - minimum number of reads to support a splice junction without the donor exon boundary conserved
	Int? variant_quality_threshold # quality threshold for variants
	#? no_missed_cleavage 		   # defaults to allow for one missed cleavage; used for generating tryptic peptides which is not yet ready

	Int? memory
  	Int? disk_space
  	Int? num_threads
  	Int? num_preemptions


	command {
		set -euo pipefail
		som_call=''
		germ_call=''
		junc_call=''
		fus_call=''

        if [ "${sep='' input_somatic_vcfs}" != '' ]; then
        	somatic="somatic"
        	mkdir $somatic
        	mv ${sep=' ' input_somatic_vcfs} $somatic
        	som_call=" --somatic somatic"
        fi
        if [ "${sep='' input_germline_vcfs}" != '' ]; then
        	germline="germline"
        	mkdir $germline
        	mv ${sep=' ' input_germline_vcfs} $germline
        	germ_call=" --germline germline"
        fi
        if [ "${sep='' input_splice_junctions_files}" != '' ]; then
        	junction="junction"
        	mkdir $junction
        	mv ${sep=' ' input_splice_junctions_files} $junction
        	junc_call=" --junction junction"
        fi
        if [ "${sep='' input_gene_fusions_files}" != '' ]; then
        	fusion="fusion"
        	mkdir $fusion
        	mv ${sep=' ' input_gene_fusions_files} $fusion
        	fus_call=" --fusion fusion"
        fi
        quilts_call="$som_call$germ_call$junc_call$fus_call"
        
        /src/QUILTS/src/quilts_references.sh \
        	-g ${reference_genome} \
        	-p ${reference_proteome}
        
		python2 /src/QUILTS/pyQUILTS/quilts.py \
		--output_dir ${output_dir} \
		--genome ./genome \
		--proteome ./proteome \
		$quilts_call \
		${"--junction_file_type " + junction_file_type} \
		${"--threshB " + threshB} \
		${"--threshD " + threshD} \
		${"--threshN " + threshN} \
		${"--variant-quality-threshold " + variant_quality_threshold}

		home=`pwd`
		cd ${output_dir}
		base=$(basename $PWD)
        cd ..
		tar -czf ${out_file} $base
        mv ${out_file} $home
	}

	output {
		File outputs = "${out_file}"
	}

	runtime {
		docker : "broadcptac/panoply_quilts:1_1"
		memory : select_first ([memory, 12]) + "GB"
		disks : "local-disk " + select_first ([disk_space, 20]) + " SSD"
		cpu : select_first ([num_threads, 1]) + ""
		preemptible : select_first ([num_preemptions, 0])
	}

	meta {
    author : "Myranda Maynard"
    email : "proteogenomics@broadinstitute.org"
  }
}

workflow quilts_workflow {
	call quilts
}