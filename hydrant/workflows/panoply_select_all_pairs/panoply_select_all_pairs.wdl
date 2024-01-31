#
# Copyright (c) 2023 The Broad Institute, Inc. All rights reserved.
#


workflow panoply_select_all_pairs {
  Array[Pair[String?,File?]]+ pairs_input

  ## Separate Array[Pair[String,File?]] into array of extant labels and extant files
  scatter (pair in pairs_input) {             ## for each pair
    if (defined(pair.left) && defined(pair.right)) {      ## if label AND file are defined
      String pair_string_ = "${pair.left}"          ## force WDL to interpret String as extant
      File pair_file_ = "${pair.right}"           ## force WDL to interpret File as extant
    }
  }

  output {
    Array[String] pair_string = select_all(pair_string_)            ## select extant labels
    Array[File] pair_file = select_all(pair_file_)                  ## select extant files

    Array[Pair[String,File]] pairs = zip(pair_string, pair_file)    ## zip extant inputs into fully-defined array
  }
}
