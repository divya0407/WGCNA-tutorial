# THIS ONLY WORKS FOR DATASETS 4107 AND 33113 - NOT ANY OTHER!!

library(hgu133plus2.db)
library(AnnotationDbi)
library(dplyr)
library(tibble)

normalized_file <- "PATH_TO_NORM_CSV"
output_file <- "PATH_TO_OUTPUT_FOLDER" # NAME THE FILE HERE ITSELF


if (!file.exists(normalized_file)) stop("Normalized file missing!")
expr_matrix <- read.csv(normalized_file, row.names = 1, check.names = FALSE)


probe_ids <- rownames(expr_matrix)
annot <- AnnotationDbi::select(hgu133plus2.db, 
                               keys = probe_ids, 
                               columns = "SYMBOL", 
                               keytype = "PROBEID")


expr_annotated <- expr_matrix %>%
  tibble::rownames_to_column("PROBEID") %>%
  dplyr::left_join(annot, by = "PROBEID") %>%
  dplyr::filter(!is.na(SYMBOL) & SYMBOL != "") %>%
  dplyr::select(-PROBEID) %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::summarise(across(everything(), mean, na.rm = TRUE)) %>%
  dplyr::ungroup()


if (exists("expr_annotated") && nrow(expr_annotated) > 0) {
  if (!dir.exists(dirname(output_file))) dir.create(dirname(output_file), recursive = TRUE)
  write.csv(expr_annotated, output_file, row.names = FALSE)
  message("SUCCESS: File produced with ", nrow(expr_annotated), " unique genes.")
} else {
  stop("ERROR: expr_annotated was not created correctly.")
}
