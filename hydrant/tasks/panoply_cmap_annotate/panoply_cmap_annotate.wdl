#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
task panoply_cmap_annotate {
    Float? ram_gb
    Int? local_disk_gb
    Int? num_preemptions

    #**Define additional inputs here**

    command {
        set -euo pipefail

        #**Command goes here**
    }

    output {
        #** Define outputs here**
    }

    runtime {
        docker : "<namespace>/panoply_cmap_annotate:1"
        memory: "${if defined(ram_gb) then ram_gb else '2'}GB"
        disks : "local-disk ${if defined(local_disk_gb) then local_disk_gb else '10'} HDD"
        preemptible : "${if defined(num_preemptions) then num_preemptions else '0'}"
    }

    meta {
        author : "D R Mani"
        email : "proteogenomics@broadinstitute.org"
    }
}

workflow panoply_cmap_annotate {

    call panoply_cmap_annotate {
        input: #**Define call inputs for panoply_cmap_annotate here**
    }

    output {
        #**Define workflow outputs here. If defined, these will be the only
        #  outputs available in the Method Configuration**
    }
}
