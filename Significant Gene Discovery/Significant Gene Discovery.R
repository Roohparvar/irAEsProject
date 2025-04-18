# # Split marker genes by cluster and export each cluster to a separate sheet in an Excel file
###################################################################################################
library(dplyr)
library(writexl)

data <- read.csv("AllMarkers.csv")

clusters <- unique(data$cluster)

sheets <- list()

for (cl in clusters) {
  sheets[[paste0("Cluster_", cl)]] <- filter(data, cluster == cl)
}

write_xlsx(sheets, path = "clustered_output.xlsx")



# Filter marker genes in each sheet based on significance thresholds and save to a new Excel file
###################################################################################################
library(openxlsx)

input_file <- "clustered_output.xlsx"
output_file <- "filtered_output.xlsx"

sheet_names <- getSheetNames(input_file)

wb <- createWorkbook()

for (sheet_name in sheet_names) {
  data <- read.xlsx(input_file, sheet = sheet_name)
  
  filtered_data <- subset(data, 
                          p_val_adj < 0.05 & 
                            avg_log2FC > 0.1 & 
                            pct.1 > 0.1 & 
                            (pct.1 - pct.2) > 0.1)
  
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet = sheet_name, x = filtered_data)
}

saveWorkbook(wb, output_file, overwrite = TRUE)

###################################################################################################