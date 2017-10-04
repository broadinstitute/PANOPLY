

task outlier_analysis {
  File inputTarball
  String outputTarball

  command {
    set -euo pipefail
    /prog/outlier.sh ${inputTarball} ${outputTarball}
  }

  output {
    File outputs = "${outputTarball}"
  }

  runtime {
    docker : "huazhou2/outlier:v2"
  }

  meta {
    author : "Hua Zhou"
    email : "Hua.Zhou@nyumc.org"
	}

}



workflow run_outlier_analysis {
  File inputFile
  String outputName

  call outlier_analysis {input: inputTarball=inputFile, outputTarball=outputName}
}

