import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_parse_sm_table/versions/1/plain-WDL/descriptor" as parse_sm_table
import "https://api.firecloud.org/ga4gh/v1/tools/broadcptac:panoply_normalize_ms_data/versions/1/plain-WDL/descriptor" as normalize


workflow panoply_process_SM_table {

  ## inputs
  String job_identifier
  String ome_type
  File sample_annotation
  File input_ssv
  File? custom_parameters


  call parse_sm_table.panoply_parse_sm_table as parse {
    input:
      SMtable = input_ssv,
      exptDesign = sample_annotation,
      analysisDir = job_identifier,
      type = ome_type,
      params = custom_parameters
  }

  call normalize.panoply_normalize_ms_data  as norm {
    input:
      inputData = parse.outputs, 
      type = ome_type,
      standalone = "false",
      params = custom_parameters
  }

  output {
    File output_tar = norm.outputs
  }

}
