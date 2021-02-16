#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#

import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_cons_clust/versions/1/plain-WDL/descriptor" as cons_clust_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_cons_clust_report/versions/1/plain-WDL/descriptor" as cons_clust_report_wdl

workflow panoply_cons_clust_workflow {
	File inputData
	String type
	File? groupsFile
	String job_identifier
	String standalone = "true"
	File yaml = '/prot/proteomics/Projects/PGDAC/src/master-parameters.yaml'

	call cons_clust_wdl.panoply_cons_clust {
    input:
      inputData = inputData,
      type = type,
      groupsFile = groupsFile,
      standalone = standalone,
      yaml = yaml,
      analysisDir = job_identifier
  }
  
  call cons_clust_report_wdl.panoply_cons_clust_report {
    input:
      tar_file = panoply_cons_clust.outputs,
      yaml_file = yaml,
      label = job_identifier,
      type = type
  }

  output {
  	File cons_clust_report = panoply_cons_clust_report.report_out
  	File cons_clust_tar = panoply_cons_clust.outputs
  }
}