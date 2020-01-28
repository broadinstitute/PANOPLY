task pgdac_aggregate_fasta {
	
	#Inputs defined here
	Array[File] snv_fasta
	Array[File] indel_fasta
	Array[File] junc_fasta
	
	command {

# file list SNVs
python <<CODE	
snv_fasta = '${sep=',' snv_fasta}'.split(",")
with open("snv_fasta.txt", "w") as fout:
	for i in range(len(snv_fasta)):
		fout.write(snv_fasta[i] + "\n")
CODE

# file list indels
python <<CODE
indel_fasta = '${sep=',' indel_fasta}'.split(",")
with open("indel_fasta.txt", "w") as fout:
	for i in range(len(indel_fasta)):
		fout.write(indel_fasta[i] + "\n")
CODE

# file list junctions
python <<CODE
junc_fasta = '${sep=',' junc_fasta}'.split(",")
with open("junc_fasta.txt", "w") as fout:
	for i in range(len(junc_fasta)):
		fout.write(junc_fasta[i] + "\n")
CODE

	}

	output {
		#Outputs defined here
		#File test="test.txt"
	}

	runtime {
		docker : "broadcptac/pgdac_cpdb:2"
	}

	meta {
		author : "Karsten Krug"
		email : "karsten@broadinstitute.org"
	}

}

workflow pgdac_aggregate_fasta_workflow {
	call pgdac_aggregate_fasta
}

