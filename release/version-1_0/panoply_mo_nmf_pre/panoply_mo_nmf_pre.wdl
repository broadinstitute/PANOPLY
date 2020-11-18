#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
################################################
## prepare input tarball for panoply_mo_nmf
##
task panoply_mo_nmf_pre {

    String label

    Array[File?] omes

    File? rna_ome
    File? cna_ome 

    Int? memory
    Int? disk_space
    Int? num_threads
    Int? num_preemptions

    command {
        set -euo pipefail
        
        prote_ome=''
        phospho_ome=''
        acetyl_ome=''
        ubiquityl_ome=''
       
        for f in ${sep=" " omes}; do
            if [ "$f" ]; then
                echo $f "FILE"
                base="$(basename $f)"
                ome=`echo $base | cut -d '-' -f 1`
                echo $ome "OME"
            else
                ome=''
            fi

            if [ "$ome" == "proteome" ];then prote_ome=$f;fi

            if [ "$ome" == "phosphoproteome" ];then phospho_ome=$f;fi

            if [ "$ome" == "acetylome" ];then acetyl_ome=$f;fi
            
            if [ "$ome" == "ubiquitylome" ];then ubiquityl_ome=$f;fi
        done
        
        echo $prote_ome "prote_ome"
        echo $phospho_ome "phospho_ome"
        echo $acetyl_ome "acetyl_ome"
        echo $ubiquityl_ome "ubi_ome"

        if [ "$prote_ome" != '' ];then echo $prote_ome > prote_ome.txt;else echo '' > prote_ome.txt;fi
        if [ "$phospho_ome" != '' ];then echo $phospho_ome > phospho_ome.txt;else echo '' > phospho_ome.txt;fi
        if [ "$acetyl_ome" != '' ];then echo $acetyl_ome > acetyl_ome.txt;else echo '' > acetyl_ome.txt;fi
        if [ "$ubiquityl_ome" != '' ];then echo $ubiquityl_ome > ubiquityl_ome.txt;else echo '' > ubiquityl_ome.txt;fi
        
        echo `cat prote_ome.txt`
        echo `cat phospho_ome.txt`
        echo `cat acetyl_ome.txt`
        echo `cat ubiquityl_ome.txt`

        echo  "${if defined(rna_ome) then rna_ome else ''}" > rna_ome.txt
        echo  "${if defined(cna_ome) then cna_ome else ''}" > cna_ome.txt

        /home/pgdac/src/create-tar.R ${label}

    }

    output {
        # Define outputs here
       File tar="${label}.tar"
    
    }

    runtime {
        docker : "broadcptac/panoply_mo_nmf_pre:1_0"
        memory: select_first ([memory, 2]) + "GB"
        disks : "local-disk " + select_first ([disk_space, 10]) + " SSD"
        cpu   : select_first ([num_threads, 1]) + ""
        preemptible : select_first ([num_preemptions, 0])
    }

    meta {
        author : "Karsten Krug"
        email : "proteogenomics@broadinstitute.org"
    }
}
## workflow
workflow panoply_mo_nmf_pre_wf {

    #File cfg_yml

    call panoply_mo_nmf_pre #{
     #   input:
     #      cfg_yml=cfg_yml
    #}

    output {
        File nmf_tar=panoply_mo_nmf_pre.tar
    }
}

