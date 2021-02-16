#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_immune_analysis/versions/21/plain-WDL/descriptor" as immune_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_immune_analysis_report/versions/2/plain-WDL/descriptor" as immune_report_wdl

workflow panoply_immune_analysis_workflow {
    File inputData
  	String type
  	File? groupsFile
    String job_identifier
    String standalone = "true"
    File yaml = '/prot/proteomics/Projects/PGDAC/src/master-parameters.yaml'

  	Float? fdr
  	Int? heatmapWidth
  	Int? heatmapHeight
    
        
    call immune_wdl.panoply_immune_analysis as immune {
    	input:
        	inputData = inputData,
          type = type,
          standalone = standalone,
          yaml = yaml,
          analysisDir = job_identifier,
          groupsFile = groupsFile,
          fdr = fdr,
          heatmapWidth = heatmapWidth,
          heatmapHeight = heatmapHeight
    }
     
	call immune_report_wdl.panoply_immune_analysis_report as immune_report {
    	input:
        	tar_file = immune.outputs,
          yaml_file = immune.yaml_file,
          label = job_identifier
    }
    
    output {
      File outputs = immune.outputs
      File report = immune_report.report_out
    }
}
