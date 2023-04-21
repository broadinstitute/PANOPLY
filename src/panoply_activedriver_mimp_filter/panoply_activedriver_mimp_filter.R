library(pacman)
p_load(dplyr)
p_load(unheadr)
p_load(stringr)

tar_files <- list.files(path = ".", pattern = "tar")
for (tar in tar_files) {
  untar(tar)
}

activedriver_results <- readRDS("activedriver_results_dir/results_all.rds")
mimp_rewiring_events <- read.csv("mimp_results_dir/mimp_output_kinase_rewiring_events_all.csv") # PART 1
mimp_gain_events <- read.csv("mimp_results_dir/all_STY_gain_mutation.csv") # PART 2
mimp_loss_events <- read.csv("mimp_results_dir/all_STY_loss_mutation.csv") # PART 2

dir.create("mimp_results_filtered_dir")

### PART 1: FILTER REWIRING EVENTS
sign_mut <- activedriver_results$merged_report %>% filter(active_region_p_adj <= 0.05)
sign_mut_mimp_filt <- sign_mut[, c("gene", "sample_id")]
colnames(sign_mut_mimp_filt) <- c("protein_id", "sample_id")
sign_mut_mimp_filt$mutation <- paste0(sign_mut$wt_residue, sign_mut$mut_position, sign_mut$mut_residue)

sign_mut_mimp_filt$uniq_mut_id <- paste0(sign_mut_mimp_filt$protein_id, "_", sign_mut_mimp_filt$mutation)
mimp_rewiring_events$uniq_mut_id <- paste0(mimp_rewiring_events$protein_id, "_", mimp_rewiring_events$mutation)
mimp_rewiring_events_filt <- mimp_rewiring_events %>% filter(uniq_mut_id %in% sign_mut_mimp_filt$uniq_mut_id)
write.csv(mimp_rewiring_events_filt, "mimp_results_filtered_dir/mimp_output_kinase_rewiring_events_activedriver_filtered.csv")


mimp_rewiring_events_filt_compact <- as.data.frame(mimp_rewiring_events_filt %>% unwrap_cols(groupingVar = uniq_mut_id, separator = "|"))

# since we unwrap based on unique combination of protein ID and mutation location, each row will only contain 1 protein ID and 1 mutation
mimp_rewiring_events_filt_compact$protein_id <- str_extract(mimp_rewiring_events_filt_compact$uniq_mut_id, "ENSP\\d+.\\d+")  # extract 1 protein ID
mimp_rewiring_events_filt_compact$mutation <- str_extract(mimp_rewiring_events_filt_compact$uniq_mut_id, "[A-Z]\\d+[A-Z]")  # extract 1 mutation

# keep most important columns
mimp_rewiring_events_filt_compact <- subset(mimp_rewiring_events_filt_compact, select = c("protein_id", "mutation", "kinase_pwm", "effect", "sample_id"))
mimp_rewiring_events_filt_compact$sample_id <-
  unlist(
    lapply( 
      lapply( # extract only unique samples
        str_split(mimp_rewiring_events_filt_compact$sample_id, pattern = "\\|"),
        unique
      ),
      paste, # stitch back with separator
      collapse = "|"
    )  
  ) # unlist back to re-assign to original sample column

write.csv(mimp_rewiring_events_filt_compact, "mimp_results_filtered_dir/mimp_output_kinase_rewiring_events_activedriver_filtered_compact.csv")
### PART 1 END


### PART 2 FILTER GAIN & LOSS EVENTS
mimp_gain_events$uniq_mut_id <- paste0(mimp_gain_events$protein_id, "_", mimp_gain_events$mutation)
mimp_loss_events$uniq_mut_id <- paste0(mimp_loss_events$protein_id, "_", mimp_loss_events$mutation)

mimp_gain_events_filt <- mimp_gain_events %>% filter(uniq_mut_id %in% sign_mut_mimp_filt$uniq_mut_id)
mimp_loss_events_filt <- mimp_loss_events %>% filter(uniq_mut_id %in% sign_mut_mimp_filt$uniq_mut_id)

write.csv(mimp_gain_events_filt, "mimp_results_filtered_dir/STY_gain_mutation_activedriver_filtered.csv")
write.csv(mimp_loss_events_filt, "mimp_results_filtered_dir/STY_loss_mutation_activedriver_filtered.csv")

### PART 2 END

