#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
workflow panoply_make_pairs_workflow {
  Array[File?] files
  String? suffix 

  scatter (m in files) {
    if ("${m}" != "") {
      String out_name = basename("${m}", "${suffix}" )
      String out_file = "${m}"
    }
  }

  output {
    Array[Pair[String?, File?]] zipped = zip(out_name, out_file)
  }
}