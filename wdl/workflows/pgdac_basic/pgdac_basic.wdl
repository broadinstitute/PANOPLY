task parse_sm_table {
	File SMtable
	File exptDesign
	String analysisDir
	String codeDir = "/prot/proteomics/Projects/PGDAC/src"
	String dataDir = "/prot/proteomics/Projects/PGDAC/data"

	command {
		set -euo pipefail
		/prot/proteomics/Projects/PGDAC/src/run-pipeline.sh inputSM -s ${SMtable} -r ${analysisDir} -c ${codeDir} -d ${dataDir} -e ${exptDesign}
	}

	output {
		File outputs = "inputSM-output.tar"
	}

	runtime {
		docker : "broadcptac/pgdac_basic:1"
	}

	meta {
		author : "D. R. Mani"
		email : "manidr@broadinstitute.org"
	}

}



task mrna_protein_corr {
	File tarball
	File rnaExpr
	String codeDir = "/prot/proteomics/Projects/PGDAC/src"

	command {
		set -euo pipefail
		/prot/proteomics/Projects/PGDAC/src/run-pipeline.sh RNAcorr -i ${tarball} -c ${codeDir} -rna ${rnaExpr}
	}

	output {
		File outputs = "RNAcorr-output.tar"
	}

	runtime {
		docker : "broadcptac/pgdac_basic:1"
	}

	meta {
		author : "D. R. Mani"
		email : "manidr@broadinstitute.org"
	}

}



workflow pgdac_basic {
	File SMtable
	File exptDesign
	File rnaExpr
	String analysisDir

  call parse_sm_table {
    input:
      SMtable=SMtable, 
      analysisDir=analysisDir, 
      exptDesign=exptDesign
  }
  
  call mrna_protein_corr {
    input: 
      tarball=parse_sm_table.outputs, 
      rnaExpr=rnaExpr
  }
  
  output {File output=mrna_protein_corr.outputs}
}

