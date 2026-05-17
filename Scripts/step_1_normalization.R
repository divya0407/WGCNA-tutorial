if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

#oligo install
if (!requireNamespace("oligo", quietly = TRUE)) {
    BiocManager::install("oligo", ask = FALSE, update = FALSE)
}

library(oligo)

# dir paths
data_dir <- "PATH_TO_DATASET"
output_dir <- "PATH_TO_OUTPUT_FOLDER"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# normalize
celFiles <- list.files(path = data_dir, full.names = TRUE, pattern = "\\.CEL\\.gz$")

if (length(celFiles) == 0) {
    stop("No CEL.gz files found! Check your path.")
}

affyRaw <- read.celfiles(celFiles)
eset <- rma(affyRaw)

# export it
expr_matrix <- exprs(eset)
output_file <- file.path(output_dir, "rma_norm_GSE113513.csv")
write.csv(expr_matrix, file = output_file)

message("Success! File saved to: ", output_file)