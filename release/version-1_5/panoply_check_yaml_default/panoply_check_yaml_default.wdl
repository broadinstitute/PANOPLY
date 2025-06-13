#
# Copyright (c) 2023 The Broad Institute, Inc. All rights reserved.
#


# Compare Terra toggle to YAML toggle, to determine if 
task panoply_check_yaml_default {
  String? param          # current (Terra set) parameter value, may or may not be set
  File yaml               # yaml file with defaults
  String param_lookup     # parameter lookup value in yaml

  command <<<
    R -s -e "if ('${param}'=='') {cat(yaml::read_yaml('${yaml}')[['${param_lookup}']])} else {cat('${param}')}"
  >>>

  output {
    Boolean param_boolean = read_boolean( stdout () )
  }
  
  

  runtime {
    docker : "broadcptacdev/panoply_common:latest"
  }
  

  meta {
    author : "C.M. Williams"
    email : "proteogenomics@broadinstitute.org"
  }
}

workflow panoply_check_yaml_default_workflow {
  call panoply_check_yaml_default
}
