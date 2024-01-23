#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
task panoply_nmf_report {
    String label

    # inputs from panoply_nmf
    File nmf_results                ## Rdata file containing results of nmf()
    Int nclust                      ## best clustering assignment
    File expr_comb                  ## combined expression GCT file
    File expr_comb_nn               ## combined expression GCT file (non negative)
    File nmf_parameters             ## Rdata file containing parameters list-object (opt) from panoply_nmf

    # inputs from panoply_nmf_postprocess
    File postprocess_tarball        ## tarball containing all figures and outputs from panply_nmf_postprocess
    File postprocess_parameters     ## Rdata file containing parameters list-object (opt) from panoply_nmf_postprocess

    Int? memory
    Int? disk_space
    Int? num_threads
    Int? num_preemptions

    command {
        set -euo pipefail
        Rscript /prot/proteomics/Projects/PGDAC/src/nmf-renderRMD.R ${"-n " + nmf_results} ${"-r " + nclust} ${"-e " + expr_comb} ${"-f " + expr_comb_nn} ${"-p " + nmf_parameters} ${"-q " + postprocess_parameters} ${"-x " + label} ${"-t " + postprocess_tarball}
    }

    output {
	File report = label + "_nmf_report.html"
    }

    runtime {
        docker : "broadcptacdev/panoply_nmf_report:latest"
        memory: select_first ([memory, 16]) + "GB"
        disks : "local-disk " + select_first ([disk_space, 10]) + " SSD"
        cpu : select_first ([num_threads, 1]) + ""
        preemptible : select_first ([num_preemptions, 0])
    }

    meta {
        author : "Karsten Krug"
        email : "proteogenomics@broadinstitute.org"
    }
}

workflow panoply_nmf_report_workflow {
    call panoply_nmf_report
}
