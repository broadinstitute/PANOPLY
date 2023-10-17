workflow panoply_scion_workflow {
    call panoply_scion
}

task panoply_scion {
    String ome
    File pome_gct_file
    File mrna_gct_file
    File TF_file

    String? dir_name
    String? type
    Float? weight_threshold
    Int? num_cores
    Boolean? verbose

    Int? memory
    Int? disk_space
    Int? num_preemptions

    command {
        set -euo pipefail

        output_dir="$(pwd)/scion_output"
  	mkdir -p $output_dir
  	cd $output_dir

        Rscript /prot/proteomics/Projects/PGDAC/src/scion/SCION.R \
        --prefix ${ome} \
        --pome.gct.file ${pome_gct_file} \
        --mrna.gct.file ${mrna_gct_file} \
        --TF.file ${TF_file} \
        --dir.name ${default=exp dir_name} \
        --type ${default=SM type} \
        --weightthreshold ${default=0 weight_threshold} \
        --num.cores ${default=1 num_cores} \
        --verbose + ${default=F verbose} \
        --libdir /prot/proteomics/Projects/PGDAC/src/scion

        tar -czvf "panoply_scion_output.tar" $(basename $out_dir)
    }

    output {
        File scion_tar = "panoply_scion_output.tar"
    }

    runtime {
        docker : "broadcptacdev/panoply_scion:latest"
        memory : select_first ([memory, 32]) + "GB"
        disks : "local-disk " + select_first ([disk_space, 50]) + " SSD"
        cpu : select_first ([num_cores, 1]) + ""
        preemptible : select_first ([num_preemptions, 0])
    }

    meta {
        author : "Natalie Clark"
        email : "nclark@broadinstitute.org"
    }
}
