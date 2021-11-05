MAIN_WD <- getwd()
project <- Sys.getenv('WORKSPACE_NAMESPACE')
workspace <- Sys.getenv('WORKSPACE_NAME')
bucket <- Sys.getenv('WORKSPACE_BUCKET')

init_project_dir <- function() {
    data_dir <<- file.path(MAIN_WD, "data")
    project_input <<- file.path(data_dir, PROJECT_NAME, "input")
    project_output <<- file.path(data_dir, PROJECT_NAME, "output")

    dir.create(data_dir, showWarnings = FALSE)
    dir.create(file.path(data_dir, PROJECT_NAME), showWarnings = FALSE)

    dir.create(project_input, showWarnings = FALSE)
    dir.create(project_output, showWarnings = FALSE)
}

list_files_in_bucket <- function(only_gct = FALSE) {
    file_paths <- system(paste0("gsutil ls ", bucket), intern=TRUE)
    files <- gsub(paste0(bucket, "/"), "", file_paths)
    if (only_gct) {
        files <- files[grepl(".gct", files)]
    }

    print(files)
}

copy_from_bucket_to_project_dir <- function(filename) {
    system(paste("gsutil cp", file.path(bucket, input_file), project_input), intern = TRUE)
}

get_ptm_sig_db <- function(id_type_out, organism) {
    if (id_type_out == "uniprot") {
        id_opt <- "uniprot"
    } else if (id_type_out == "seqwin") {
        id_opt <- "flanking"
    } else if (id_type_out == "psp") {
        id_opt <- "sitegrpid"
    } else {
        print("unsupported `id_type_out` was selected")
    }

    ptm_sig_db <- file.path("~/db/ptmsigdb", paste0("ptm.sig.db.all.", id_opt, ".human.v1.9.0.gmt"))
    return(ptm_sig_db)
}

preprocess_gct <- function() {
    setwd(project_output)  # change to folder to write

    sys_out <- system(paste("~/src/preprocessGCT.R",
        "-i", input_ds,
        "-l", "ssc",  # single-site centric
        "-t", id_type,
        "-o", id_type_out,
        "-a", acc_type_in,
        "-s", seqwin_col,
        "-d", localized,
        "-m", mode,
        "-r", residue,
        "-p", ptm,
        "-u", TRUE,
        "-z", "~/src/"
    ), intern = TRUE)

    print(sys_out)
    setwd(MAIN_WD)  # change to normal work directory
    fn <<- file.path(project_output, "fn.out")
    input_ds_proc <<- file.path(project_output, system(paste("cat", fn), intern = TRUE))  # locate processed GCT file
}

run_ptm_sea <- function() {
    setwd(project_output)  # change to folder to write

    sys_out <- system(paste("~/ssgsea-cli.R",
        "-i", input_ds_proc,
        "-d", ptm_sig_db,
        "-o", output_prefix,
        "-n", sample_norm_type,
        "-w", weight,
        "-c", correl_type,
        "-t", statistic,
        "-s", output_score_type,
        "-p", nperm,
        "-m", min_overlap,
        "-x", extended_output,
        "-e", export_signal_gct,
        "-g", global_fdr,
        "-l", TRUE
    ), intern = TRUE)

    print(sys_out)
    setwd(MAIN_WD)  # change to normal work directory
}

save_results_to_bucket <- function() {
    setwd(project_output) 
    setwd("../")

    output_zip <- paste0(PROJECT_NAME, ".zip")
    sys_out <- system(paste("zip -r", output_zip, "output"), intern = TRUE)
    system(paste("gsutil cp", output_zip, file.path(bucket, output_zip)), intern = TRUE)

    print(sys_out)
    setwd(MAIN_WD)  # change to normal work directory
}
