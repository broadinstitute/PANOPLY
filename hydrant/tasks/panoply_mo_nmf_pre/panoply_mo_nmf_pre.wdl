################################################
## prepare input tarball for panoply_mo_nmf
##
task panoply_mo_nmf_pre {

    String label

    File? prote_ome
    File? phospho_ome
    File? acetyl_ome
    File? ubiquityl_ome

    File? rna_ome
    File? cna_ome 

    Float? ram_gb
    Int? local_disk_gb
    Int? num_preemptions

    command {
        set -euo pipefail
	
        echo  "${if defined(prote_ome) then prote_ome else ''}" > prote_ome.txt
	echo  "${if defined(phospho_ome) then phospho_ome else ''}" > phospho_ome.txt
 	echo  "${if defined(acetyl_ome) then acetyl_ome else ''}" > acetyl_ome.txt
 	echo  "${if defined(ubiquityl_ome) then ubiquityl_ome else ''}" > ubiquityl_ome.txt

        echo  "${if defined(rna_ome) then rna_ome else ''}" > rna_ome.txt
	echo  "${if defined(cna_ome) then cna_ome else ''}" > cna_ome.txt

        /home/pgdac/src/create-tar.R ${label}

    }

    output {
        # Define outputs here
	File tar="${label}.tar"
	
    }

    runtime {
        docker : "broadcptacdev/panoply_mo_nmf_pre:latest"
        memory: "${if defined(ram_gb) then ram_gb else '2'}GB"
        disks : "local-disk ${if defined(local_disk_gb) then local_disk_gb else '10'} HDD"
        preemptible : "${if defined(num_preemptions) then num_preemptions else '0'}"
    }

    meta {
        author : "Karsten Krug"
        email : "karsten@broadinstitute.org"
    }
}

## workflow
workflow panoply_mo_nmf_pre_workflow {

    String label

    call panoply_mo_nmf_pre #{
        input:
     	    label=label
    }

    output {
        File nmf_tar=panoply_mo_nmf_pre.tar
    }
}
