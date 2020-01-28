

task outlier_analysis {
  File inputTarball
  String outputTarball
  Int plotTopN

  command {
    set -euo pipefail
    /prog/outlier.sh ${inputTarball} ${outputTarball} ${plotTopN}
  }

  output {
    File outputs = "${outputTarball}"
  }

  runtime {
    docker : "huazhou2/outlier:v4"
  }

  meta {
    author : "Hua Zhou"
    email : "Hua.Zhou@nyumc.org"
  }

}



workflow run_outlier_analysis {
  File inputFile
  String outputName
  Int plotTopN

  call outlier_analysis {input: inputTarball=inputFile, outputTarball=outputName, plotTopN=plotTopN}
}

