#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_blacksheep/versions/11/plain-WDL/descriptor" as blacksheep_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_blacksheep_report/versions/6/plain-WDL/descriptor" as blacksheep_report_wdl

workflow panoply_blacksheep_workflow {
    File input_gct
    File master_yaml
    String output_prefix
    File? groups_file
    String type
    
    call blacksheep_wdl.panoply_blacksheep {
        input:
            input_gct = input_gct,
            master_yaml = master_yaml,
            output_prefix = output_prefix,
            groups_file = groups_file
    }
    
    call blacksheep_report_wdl.panoply_blacksheep_report {
        input:
            input_tar = panoply_blacksheep.tar_out,
            output_prefix = output_prefix,
            type = type
    }
    
    output{
        File blacksheep_tar = panoply_blacksheep.tar_out
        File blacksheep_report = panoply_blacksheep_report.report_out
    }
        
}
