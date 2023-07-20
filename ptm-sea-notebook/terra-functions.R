script.dir <<- "/ptm-sea"  # needed for ssGSEA2.0

source("/ptm-sea/src/ssGSEA2.0.R")
source("/ptm-sea/src/panoply_ptmsea_functions.R")


MAIN_WD <- getwd()
project <- Sys.getenv('WORKSPACE_NAMESPACE')
workspace <- Sys.getenv('WORKSPACE_NAME')
bucket <- Sys.getenv('WORKSPACE_BUCKET')


#' @title Set up directories
#'
#' @param name (optional): supply a custom PTM-SEA run name, otherwise dates
init_project_dir <- function(name = NULL) {
    data_dir <<- file.path(MAIN_WD, "data")
    project_name <<- get_run_name(name)
    project_dir <<- file.path(data_dir, project_name)
    project_input <<- file.path(project_dir, "input")
    project_output <<- file.path(project_dir, "output")

    dir.create(data_dir, showWarnings = FALSE)
    dir.create(file.path(data_dir, project_name), showWarnings = FALSE)

    dir.create(project_input, showWarnings = FALSE)
    dir.create(project_output, showWarnings = FALSE)
}

get_run_name <- function(name = NULL) {
    if (!is.null(name)) {
        run_name <- name
    } else {
        run_name <- format(Sys.time(), "%y%m%d-%H%M%S", tz = "EST")
    }
    
    return(run_name)
}

list_files_in_bucket <- function(only_gct = FALSE, only_gmt = FALSE) {
    file_paths <- system(paste0("gsutil ls -r ", bucket), intern=TRUE)
    files <- gsub(paste0(bucket, "/"), "", file_paths)
    if (only_gct) {
        files <- files[grepl(".gct", files)]
    }
    if (only_gmt) {
        files <- files[grepl(".gmt", files)]
    }

    print(files)
}

#' @title Use gsutil to copy from bucket to current environment
#'
#' @param filename: file name in the bucket (including folder path)
copy_from_bucket_to_project_dir <- function(filename) {
    system(paste("gsutil cp", file.path(bucket, filename), project_input), intern = TRUE)
}

#' @title Get correct PTMsigDB for PTM-SEA
#'
#' @param id_type_out: type of IDs in PTM (options: "uniprot", "flanking", "sitegrpid")
#' @param organism: organism database (options: "human", "mouse", "rat")
#'
#' @return protein | transcript | gene_id |
get_ptm_sig_db <- function(id_type_out, organism) {
    if (is.null(ptm_sig_db_path)) {
        if (id_type_out == "uniprot") {
            id_opt <- "uniprot"
        } else if (id_type_out == "seqwin") {
            id_opt <- "flanking"
        } else if (id_type_out == "psp") {
            id_opt <- "sitegrpid"
        } else {
            print("unsupported `id_type_out` was selected")
        }
        ptm_sig_db <- file.path("/ptm-sea/db/ptmsigdb", paste0("ptm.sig.db.all.", id_opt, ".", organism, ".v2.0.0.gmt"))
    } else {
        copy_from_bucket_to_project_dir(ptm_sig_db_path)
        ptm_sig_db <- file.path(project_input, basename(ptm_sig_db_path))
    }

    return(ptm_sig_db)
}

#' @title Runs preprocess-GCT.R to keep fully localized sites
#'
#' @return input_ds_proc: processed GCT input for PTM-SEA
preprocess_gct <- function() {
    setwd(project_output)  # change to folder to write

    out <- preprocessGCT(
        gct.str = input_ds,
        level = "ssc",  # single-site centric
        id.type = id_type,
        id.type.out = id_type_out,
        acc.type = acc_type_in,
        seqwin.col = seqwin_col,
        gene.col = gene_symbol_col,
        humanize.gene.names = humanize_gene,
        loc = localized,
        mode = mode,
        mod.res = residue,
        mod.type = ptm,
        appenddim = FALSE,
        preprocess.gct = TRUE
    )

    out_path <- file.path(project_output, "fn.out")
    writeLines(out, con = out_path)
    input_ds_proc <<- file.path(project_output, system(paste("cat", out_path), intern = TRUE))  # locate processed GCT file

    setwd(MAIN_WD)  # change to normal work directory
}

#' @title Runs PTM-SEA on preprocessed 
#'
#' @return input_ds_proc
run_ptm_sea <- function(save_to_bucket = TRUE, name = NULL) {
    setwd(project_output)  # change to folder to write
    
    log.file <- file.path(project_output, paste0(output_prefix, '_ssgsea.log.txt'))
    res <- ssGSEA2(
        input.ds = input_ds_proc,
        gene.set.databases = ptm_sig_db,
        output.prefix = output_prefix,
        sample.norm.type = sample_norm_type,
        weight = weight,
        correl.type = correl_type,
        statistic = statistic,
        output.score.type = output_score_type,
        nperm = nperm,
        min.overlap = min_overlap,
        extended.output = extended_output,
        export.signat.gct = export_signal_gct,
        global.fdr = global_fdr,
        par = TRUE,
        spare.cores = 0,
        log.file = log.file
    )

    setwd(MAIN_WD)  # change to normal work directory
    
    if (save_to_bucket) {
        save_results_to_bucket(name)
    }
}
#' @title Zips outputs and saves back to original bucket
save_results_to_bucket <- function(name = NULL) {
    setwd(project_dir) 
    
    if (!is.null(name)) {
        output_zip <- paste0(name, ".zip")
    } else {
        output_zip <- paste0(project_name, ".zip")
    }

    sys_out <- system(paste("zip -r", output_zip, "output"), intern = TRUE, ignore.stdout = TRUE)
    print(paste("PTM-SEA outputs compressed:", file.path(project_dir, output_zip)))
    system(paste("gsutil cp", output_zip, file.path(bucket, output_zip)), intern = TRUE)
    print(paste("Output zip copied to the bucket:", output_zip))

    print(sys_out)
    setwd(MAIN_WD)  # change to normal work directory
}
