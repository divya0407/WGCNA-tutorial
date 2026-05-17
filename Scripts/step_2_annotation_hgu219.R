# ONLY WORKS FOR 44076

if (!requireNamespace("hgu219.db", quietly = TRUE)) {
  BiocManager::install("hgu219.db", ask = FALSE, update = FALSE)
}
library(hgu219.db)
library(AnnotationDbi)
library(dplyr)
library(tibble)

# -------------------------------
# 1) Paths
# -------------------------------
normalized_file <- "PATH_TO_NORM_CSV"
output_file <- "PATH_TO_OUTPUT_FOLDER/GSE44076_annotated.csv" # NAME THE FILE

# 2) Load Data
if (!file.exists(normalized_file)) stop("Normalized file missing!")
expr_matrix <- read.csv(normalized_file, row.names = 1, check.names = FALSE)

# 3) Map Probes using the hgu219 dictionary
probe_ids <- rownames(expr_matrix)
annot <- AnnotationDbi::select(hgu219.db, 
                               keys = probe_ids, 
                               columns = "SYMBOL", 
                               keytype = "PROBEID")

# 4) Merge & Collapse (WGCNA Prep)
temp_data <- tibble::rownames_to_column(expr_matrix, "PROBEID")

expr_annotated <- temp_data %>%
  dplyr::inner_join(annot, by = "PROBEID") %>%
  dplyr::filter(!is.na(SYMBOL) & SYMBOL != "") %>%
  dplyr::select(-PROBEID) %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::summarise(across(everything(), mean, na.rm = TRUE)) %>%
  dplyr::ungroup()

# 5) Save
if (!dir.exists(dirname(output_file))) dir.create(dirname(output_file), recursive = TRUE)
write.csv(expr_annotated, output_file, row.names = FALSE)

message("SUCCESS: GSE44076 (HG-U219) is annotated and collapsed.")