#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
################################################
## prepare input tarball for panoply_so_nmf
##
task panoply_so_nmf_pre {

    String label

    File ome
    String ome_type
    
    Int? memory
    Int? disk_space
    Int? num_threads
    Int? num_preemptions

    command {
        set -euo pipefail
        
        # move $ome to PWD, so the directory structure isn't retained
        cp ${ome} .
        ome_fname=$(basename ${ome})
        
        # create nmf.conf
        echo -e ${ome_type}'\t'$ome_fname > nmf.conf
        
        # create tar file with nmf.conf and .gct
        tar -cvf "${label}_${ome_type}.tar" nmf.conf $ome_fname

    }

    output {
       # Define outputs here
       File tar="${label}_${ome_type}.tar"
     
    }

    runtime {
        docker : "broadcptacdev/panoply_mo_nmf_pre:DEV"
        memory: select_first ([memory, 16]) + "GB"
        disks : "local-disk " + select_first ([disk_space, 10]) + " SSD"
        cpu   : select_first ([num_threads, 1]) + ""
        preemptible : select_first ([num_preemptions, 0])
    }

    meta {
        author : "C. Williams"
        email : "proteogenomics@broadinstitute.org"
    }
}

## workflow
workflow panoply_so_nmf_pre_wf {

    call panoply_so_nmf_pre 

    output {
        File nmf_tar=panoply_so_nmf_pre.tar
    }
}

