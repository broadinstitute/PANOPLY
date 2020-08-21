import "https://api.firecloud.org/ga4gh/v1/tools/broadcptacdev:panoply_main/versions/2/plain-WDL/descriptor" as main
# import panoply_blacksheep
# import panoply_mo_nmf


task gather_omes_results {
  # task to gather all outputs from the panoply_main scatter calls
  # keep files and/or create tar files for use downstream in the workflow 
}


workflow panoply_unified_workflow {
  File? prote_ome
  File? phospho_ome
  File? acetyl_ome
  File? ubiquityl_ome
  String job_id
  ## add optional variable for all panoply_main required arguments
  
  Array[File] omes = ["${prote_ome}", "${phospho_ome}", "${acetyl_ome}", "${ubiquityl_ome}"]
  
  scatter (p in omes) {
    if ("${p}" != "") {
      String ome = basename ("${p}", "-aggregate.gct")
      String ptmsea = if "${ome}"=="phosphoproteome" then "TRUE" else "FALSE"
        call main.panoply_main as pome {
        input:
            ## include all required arguments from above
            input_pome="${p}",
            ome_type="${ome}",
            job_identifier="${job_id}-${ome}",
            run_ptmsea="${ptmsea}"
      }
    }
  }
  
  ## call gather_omes_results 
  ##   with appropriate inputs, and set up proper outputs
  ##   Array[File?] pome.pgdac_full (tar file) will be one of the inputs
  
  
  ## run panoply_blacksheep (imported as a WDL) on all input data gct files using a scatter
  ## gather blacksheep output (using a new task, if needed)
  
  
  ## assemble data tables and/or tar file for panoply_mo_nmf module
  ## call panoply_mo_nmf (imported as a WDL)

  ## assemble final output combining results from panoply_main, _blacksheep and _mo_nmf
  
  output {
    ## final output file(s)
  }
}
