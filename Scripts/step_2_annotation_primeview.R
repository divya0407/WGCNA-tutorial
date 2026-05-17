# ONLY WORKS FOR 113513

library(biomaRt)
library(tidyverse)

normalized_file <- "PATH_TO_NORM_FILE"
expr_matrix <- read.csv(normalized_file, row.names = 1, check.names = FALSE)
probe_ids <- rownames(expr_matrix)

# 2. Connect to Ensembl (Affymetrix PrimeView Map)
# This step requires internet connection
cat("Connecting to Ensembl to map PrimeView probes...\n")
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")


annot_lookup <- getBM(attributes = c("affy_primeview", "external_gene_name"),
                      filters = "affy_primeview",
                      values = probe_ids,
                      mart = mart)


colnames(annot_lookup) <- c("PROBEID", "SYMBOL")


cat("Merging and collapsing duplicates...\n")
temp_data <- expr_matrix %>% rownames_to_column("PROBEID")

expr_annotated <- temp_data %>%
  inner_join(annot_lookup, by = "PROBEID") %>%
  filter(!is.na(SYMBOL) & SYMBOL != "") %>%
  select(-PROBEID) %>%
  group_by(SYMBOL) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  ungroup()


output_file <- "PATH_TO_OUTPUT_FOLDER"
if (!dir.exists(dirname(output_file))) dir.create(dirname(output_file), recursive = TRUE)
write.csv(expr_annotated, output_file, row.names = FALSE)

cat("SUCCESS: GSE113513 is now annotated and saved.")
print(paste("Final Gene Count:", nrow(expr_annotated)))