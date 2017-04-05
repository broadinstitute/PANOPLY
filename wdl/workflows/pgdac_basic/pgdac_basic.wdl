task initialize {
	File SMtable
	File exptDesign
	String analysisDir
	String codeDir = "/prot/proteomics/Projects/PGDAC/src"

	command {
		set -euo pipefail
		/prot/proteomics/Projects/PGDAC/src/run-pipeline.sh init -s ${SMtable} -r ${analysisDir} -c ${codeDir} -e ${exptDesign}
	}

	output {
		File outputs = "init-output.tar"
	}

	runtime {
		docker : "broadcptac/pgdac_basic:1"
	}

	meta {
		author : "D. R. Mani"
		email : "manidr@broadinstitute.org"
	}

}



task parse_sm_table {
	File tarball

	command {
		set -euo pipefail
		/prot/proteomics/Projects/PGDAC/src/run-pipeline.sh parseSM -i ${tarball}
	}

	output {
		File outputs = "parseSM-output.tar"
	}

	runtime {
		docker : "broadcptac/pgdac_basic:1"
	}

	meta {
		author : "D. R. Mani"
		email : "manidr@broadinstitute.org"
	}

}



task normalize_data {
	File tarball

	command {
		set -euo pipefail
		/prot/proteomics/Projects/PGDAC/src/run-pipeline.sh norm -i ${tarball}
	}

	output {
		File outputs = "norm-output.tar"
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
	String analysisDir

  call initialize {input:SMtable=SMtable, analysisDir=analysisDir, exptDesign=exptDesign}
  call parse_sm_table {input: tarball=initialize.outputs}
  call normalize_data {input: tarball=parse_sm_table.outputs}
}

