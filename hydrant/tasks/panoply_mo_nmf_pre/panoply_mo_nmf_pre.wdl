#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
################################################
## prepare input tarball for panoply_mo_nmf
##
task panoply_mo_nmf_pre {

    String label

    Array[File?] omes
    
    Boolean? balance_omes # default is TRUE
    Float? tol  # Tolerance specifying the maximal accepted difference (as a fraction of total variance) between contributions from different data types. Used as stopping criterion to end optimization.
    Float? var  # Explained variance by PCA. Used to extract the number of PCs explaining the specified fraction of variance in the multiomics data matrix.
    String? zscore_mode # z-score mode: row (z-score rows), col (z-score columns), rowcol (z-score rows and then columns) 

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
        glyco_ome=''
       
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
            
            if [ "$ome" == "glycoproteome" ];then glyco_ome=$f;fi
            
        done
        
        echo $prote_ome "prote_ome"
        echo $phospho_ome "phospho_ome"
        echo $acetyl_ome "acetyl_ome"
        echo $ubiquityl_ome "ubi_ome"
        echo $glyco_ome "glyco_ome"

        if [ "$prote_ome" != '' ];then echo $prote_ome > prote_ome.txt;else echo '' > prote_ome.txt;fi
        if [ "$phospho_ome" != '' ];then echo $phospho_ome > phospho_ome.txt;else echo '' > phospho_ome.txt;fi
        if [ "$acetyl_ome" != '' ];then echo $acetyl_ome > acetyl_ome.txt;else echo '' > acetyl_ome.txt;fi
        if [ "$ubiquityl_ome" != '' ];then echo $ubiquityl_ome > ubiquityl_ome.txt;else echo '' > ubiquityl_ome.txt;fi
        if [ "$glyco_ome" != '' ];then echo $glyco_ome > glyco_ome.txt;else echo '' > glyco_ome.txt;fi
        
        echo `cat prote_ome.txt`
        echo `cat phospho_ome.txt`
        echo `cat acetyl_ome.txt`
        echo `cat ubiquityl_ome.txt`
        echo `cat glyco_ome.txt`

        echo  "${if defined(rna_ome) then rna_ome else ''}" > rna_ome.txt
        echo  "${if defined(cna_ome) then cna_ome else ''}" > cna_ome.txt


        ## run balanace filter
        if [ ${default="true" balance_omes} = true ];then
        
                /home/pgdac/src/filter-gcts-to-balance-omes.R -d . -t ${default="0.01" tol} -v ${default="0.9" var} -z ${default="rowcol" zscore_mode}
               
                ## update file paths
                if [ "$prote_ome" != '' ];then
                    Rscript -e "ome <-'prote_ome.txt'; fn <- readLines(ome); fn <- sub('.*/', '', fn); fn <- sub('.gct$', '-balanced-contrib.gct', fn); writeLines(fn, con=ome)"
                fi
                if [ "$phospho_ome" != '' ];then
                    Rscript -e "ome <-'phospho_ome.txt'; fn <- readLines(ome); fn <- sub('.*/', '', fn); fn <- sub('.gct$', '-balanced-contrib.gct', fn); writeLines(fn, con=ome)"
                fi
                
                if [ "$acetyl_ome" != '' ];then
                    Rscript -e "ome <-'acetyl_ome.txt'; fn <- readLines(ome); fn <- sub('.*/', '', fn); fn <- sub('.gct$', '-balanced-contrib.gct', fn); writeLines(fn, con=ome)"
                fi
                
                if [ "$ubiquityl_ome" != '' ];then
                    Rscript -e "ome <-'ubiquityl_ome.txt'; fn <- readLines(ome); fn <- sub('.*/', '', fn); fn <- sub('.gct$', '-balanced-contrib.gct', fn); writeLines(fn, con=ome)"
                fi
                
                if [ "$glyco_ome" != '' ];then
                    Rscript -e "ome <-'glyco_ome.txt'; fn <- readLines(ome); fn <- sub('.*/', '', fn); fn <- sub('.gct$', '-balanced-contrib.gct', fn);  writeLines(fn, con=ome)"
                fi
                ## cna & rna
                Rscript -e "ome <-'rna_ome.txt'; fn <- readLines(ome); fn <- sub('.*/', '', fn); fn <- sub('.gct$', '-balanced-contrib.gct', fn);  writeLines(fn, con=ome)"
                Rscript -e "ome <-'cna_ome.txt'; fn <- readLines(ome); fn <- sub('.*/', '', fn); fn <- sub('.gct$', '-balanced-contrib.gct', fn);  writeLines(fn, con=ome)"
            else
                Rscript -e "pdf('balance-omes.pdf', 10, 5);plot.new();legend('top', legend=c('Balance filter was not run!'), bty='n', cex=2);dev.off()"
        fi
                
        ##  create tar
        /home/pgdac/src/create-tar.R ${label}

    }

    output {
        # Define outputs here
       File tar="${label}.tar"
       File pdf="balance-omes.pdf"
     
    }

    runtime {
        docker : "broadcptacdev/panoply_mo_nmf_pre:latest"
        memory: select_first ([memory, 16]) + "GB"
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

    call panoply_mo_nmf_pre 

    output {
        File nmf_tar=panoply_mo_nmf_pre.tar
        File balance_omes_pdf=panoply_mo_nmf_pre.pdf
    }
}

