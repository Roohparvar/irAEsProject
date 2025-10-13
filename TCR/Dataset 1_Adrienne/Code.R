library(scRepertoire)

base_path <- "C:/Esmaeil/irAEsProject/Backup/Part 1/Data/TCR_Data/All"

subfolders <- list.dirs(base_path, recursive = FALSE, full.names = TRUE)

TCR_data <- list()

for (folder in subfolders) {
  file_path <- file.path(folder, "all_contig_annotations.csv")
  
  if (file.exists(file_path)) {
    folder_name <- basename(folder) 
    TCR_data[[folder_name]] <- read.csv(file_path, stringsAsFactors = FALSE)
  } else {
    message("File not found: ", file_path)
  }
}


TCR_data <- loadContigs(input = TCR_data)
paired_TCR <- combineTCR(TCR_data, samples = names(TCR_data), filterMulti = TRUE)

combined_df <- do.call(rbind, paired_TCR)
combined_matrix <- as.matrix(combined_df)


save_path <- "C:/Esmaeil/irAEsProject/irAEsProject/TCR/Dataset 1_Adrienne/combined_metadata_matrix.rds"

# Save combined matrix
saveRDS(combined_matrix, file = save_path)