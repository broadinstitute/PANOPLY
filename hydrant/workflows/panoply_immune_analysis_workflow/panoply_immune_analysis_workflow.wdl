#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
version 1.1

import "../../tasks/panoply_immune_analysis/panoply_immune_analysis.wdl" as immune_wdl
import "../../tasks/panoply_immune_analysis_report/panoply_immune_analysis_report.wdl" as immune_report_wdl

workflow panoply_immune_analysis_workflow {
    input {
      File inputData
    	String type
    	String standalone
    	File yaml
    	String? analysisDir
    	File? groupsFile
    	Float? fdr
    	Int? heatmapWidth
    	Int? heatmapHeight
      String label
    }

    call immune_wdl.panoply_immune_analysis as immune {
    	input:
        	inputData = inputData,
            type = type,
            standalone = standalone,
            yaml = yaml,
            analysisDir = analysisDir,
            groupsFile = groupsFile,
            fdr = fdr,
            heatmapWidth = heatmapWidth,
            heatmapHeight = heatmapHeight
    }
     
	call immune_report_wdl.panoply_immune_analysis_report as immune_report {
    	input:
        	tar_file = immune.outputs,
            yaml_file = immune.yaml_file,
            label = label
    }
    
    output {
    	File outputs = immune.outputs
        File report = immune_report.report_out
    }
}
