task panoply_blacksheep {
    Float? ram_gb
    Int? local_disk_gb
    Int? num_preemptions

    File input_gct
    String genesymbol_col
    String identifiers_file
    String annotations_cols
    Int fraction_cutoff
    Int fdr_cutoff
    String heatmap_annotations_cols

    command {
        set -euo pipefail

        Rscript blacksheep_rcode.R "${input_gct}" "${genesymbol_col}" "${identifiers_file}" "${annotations_cols}" ${fraction_cutoff} ${fdr_cutoff} "${heatmap_annotations_cols}"

        tar -czvf blacksheep_outlier_analysis.tar.gz blacksheep 
    }

    output {
        File output = "blacksheep_outlier_analysis.tar.gz"
    }

    runtime {
        docker : "broadcptacdev/panoply_blacksheep:latest"
        memory: "${if defined(ram_gb) then ram_gb else '2'}GB"
        disks : "local-disk ${if defined(local_disk_gb) then local_disk_gb else '10'} HDD"
        preemptible : "${if defined(num_preemptions) then num_preemptions else '0'}"
    }

    meta {
        author : "Karen Christianson"
        email : "karen@broadinstitute.org"
    }
}

workflow panoply_blacksheep {

    call panoply_blacksheep

}
