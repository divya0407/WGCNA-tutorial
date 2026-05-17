# =========================================================
#                INTEGRATED WGCNA PIPELINE
#          Colorectal Cancer Multi-GEO Analysis
# =========================================================

# =========================
# 1. LOAD LIBRARIES
# =========================

setwd("your/directory")

if (!require("BiocManager")) install.packages("BiocManager")

packages <- c("WGCNA", "limma", "sva", "tidyverse")

for(pkg in packages){
  if(!require(pkg, character.only = TRUE)){
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

options(stringsAsFactors = FALSE)
enableWGCNAThreads()

# =========================
# 2. LOAD DATASETS
# =========================

df1 <- read.csv("GSE4107_WGCNA_Ready.csv",
                row.names = 1,
                check.names = FALSE)

df2 <- read.csv("GSE33113_WGCNA_Ready.csv",
                row.names = 1,
                check.names = FALSE)

df3 <- read.csv("GSE44076_WGCNA_Ready.csv",
                row.names = 1,
                check.names = FALSE)

df4 <- read.csv("GSE113513_WGCNA_Ready.csv",
                row.names = 1,
                check.names = FALSE)

# =========================
# 3. FIND COMMON GENES
# =========================

common_genes <- Reduce(intersect,
                       list(rownames(df1),
                            rownames(df2),
                            rownames(df3),
                            rownames(df4)))

length(common_genes) #18370

# =========================
# 4. SUBSET COMMON GENES
# =========================

expr1_sub <- df1[common_genes, ]
expr2_sub <- df2[common_genes, ]
expr3_sub <- df3[common_genes, ]
expr4_sub <- df4[common_genes, ]

# =========================
# 5. MERGE DATASETS
# =========================

merged_expr <- cbind(expr1_sub,
                     expr2_sub,
                     expr3_sub,
                     expr4_sub)

dim(merged_expr) # 18370   392

# =========================
# 6. CREATE BATCH VECTOR
# =========================

batch <- c(rep("GSE4107", ncol(expr1_sub)),
           rep("GSE33113", ncol(expr2_sub)),
           rep("GSE44076", ncol(expr3_sub)),
           rep("GSE113513", ncol(expr4_sub)))

table(batch)

# =========================
# 7. CREATE GROUP VECTOR
# =========================

group <- c(rep("Tumor",12), rep("Normal",10),
           rep("Normal",6), rep("Tumor",90),
           rep("Normal",148), rep("Tumor",98),
           rep("Normal",14), rep("Tumor",14))

length(group) #392

# =========================
# 8. COMBAT BATCH CORRECTION
# =========================

mod <- model.matrix(~group)

merged_corrected <- ComBat(
  dat = as.matrix(merged_expr),
  batch = batch,
  mod = mod,
  par.prior = TRUE
)

# =========================
# 9. TRANSPOSE FOR WGCNA
# Rows = Samples
# Columns = Genes
# =========================

datExpr <- t(merged_corrected)

dim(datExpr) # 392 18370

# =========================
# 10. QUALITY CHECK
# =========================

gsg <- goodSamplesGenes(datExpr,
                        verbose = 3)

if (!gsg$allOK) {
  
  datExpr <- datExpr[gsg$goodSamples,
                     gsg$goodGenes]
}

# =========================
# 11. REMOVE LOW VARIANCE GENES
# =========================

geneVariance <- apply(datExpr,
                      2,
                      var)

cutoff <- quantile(geneVariance,
                   0.25)

datExpr <- datExpr[, geneVariance > cutoff]

dim(datExpr)

library(WGCNA)

options(stringsAsFactors = FALSE)

enableWGCNAThreads(nThreads = 2)

# =========================
# 12. SAMPLE CLUSTERING
# =========================

sampleTree <- hclust(dist(datExpr),
                     method = "average")

# =========================
# 13. TRAIT DATA
# =========================

traitData <- data.frame(
  Tumor  = ifelse(group == "Tumor", 1, 0),
  Normal = ifelse(group == "Normal", 1, 0)
)

rownames(traitData) <- rownames(datExpr)

# =========================
# 14. SAMPLE DENDROGRAM
# =========================

traitColors <- numbers2colors(traitData,
                              colors = c("white", "red"),
                              signed = FALSE)

pdf("SampleDendrogram_TraitHeatmap.pdf", width = 10, height = 6)

plotDendroAndColors(sampleTree,
                    traitColors,
                    groupLabels = colnames(traitData),
                    main = "Sample dendrogram and trait heatmap",
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05)

dev.off()

# =========================
# 15. SOFT THRESHOLD
# =========================

powers <- c(1:20)

sft <- pickSoftThreshold(
  datExpr,
  powerVector = powers,
  verbose = 5
)

# =========================
# 16. SCALE FREE TOPOLOGY
# =========================

pdf("2_Scale_Independence.pdf",
    width = 7,
    height = 6)

plot(
  sft$fitIndices[,1],
  -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
  xlab = "Soft Threshold (power)",
  ylab = "Scale Free Topology Model Fit, signed R^2",
  type = "n",
  main = "Scale Independence"
)

text(
  sft$fitIndices[,1],
  -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
  labels = powers,
  col = "red"
)

abline(h = 0.90,
       col = "blue")

dev.off()

# =========================
# 17. MEAN CONNECTIVITY
# =========================

pdf("3_Mean_Connectivity.pdf",
    width = 7,
    height = 6)

plot(
  sft$fitIndices[,1],
  sft$fitIndices[,5],
  xlab = "Soft Threshold (power)",
  ylab = "Mean Connectivity",
  type = "n",
  main = "Mean Connectivity"
)

text(
  sft$fitIndices[,1],
  sft$fitIndices[,5],
  labels = powers,
  col = "red"
)

dev.off()

# =========================
# 18. SELECT SOFT POWER
# =========================

softPower <- 6

# =========================
# 19. TOP VARIABLE GENES
# =========================

geneVariance <- apply(datExpr,
                      2,
                      var)

topGenes <- names(
  sort(geneVariance,
       decreasing = TRUE)
)[1:5000]

datExpr_filtered <- datExpr[, topGenes]

dim(datExpr_filtered) # 392 5000

# =========================
# 20. CONSTRUCT NETWORK
# =========================

net <- blockwiseModules(
  datExpr_filtered,
  power = softPower,
  TOMType = "signed",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  deepSplit = 2,
  numericLabels = FALSE,
  pamRespectsDendro = FALSE,
  saveTOMs = FALSE,
  verbose = 3
)

# =========================
# 21. MODULE COLORS
# =========================

moduleColors <- net$colors

table(moduleColors)

# =========================
# 22. GENE DENDROGRAM
# =========================
geneTree <- net$dendrograms[[1]]

pdf("4_Gene_Dendrogram.pdf", width = 12, height = 8)
plotDendroAndColors(geneTree,
                    moduleColors[net$blockGenes[[1]]],
                    "Module Colors",
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05)
dev.off()
# =========================
# 23. MODULE EIGENGENES
# =========================

MEs <- moduleEigengenes(
  datExpr_filtered,
  moduleColors
)$eigengenes

MEs <- orderMEs(MEs)

# =========================
# 24. MODULE TRAIT CORRELATION
# =========================

moduleTraitCor <- cor(MEs, traitData, use = "p")

moduleTraitPvalue <- corPvalueStudent(
  moduleTraitCor,
  nrow(datExpr_filtered)
)

textMatrix <- paste(
  signif(moduleTraitCor, 2),
  "\n(",
  signif(moduleTraitPvalue, 1),
  ")",
  sep = ""
)

dim(textMatrix) <- dim(moduleTraitCor)

# =========================
# 25. HEATMAP
# =========================

pdf("5_Module_Trait_Relationships.pdf",
    width = 9,
    height = 7)

par(mar = c(8, 10, 3, 3))

labeledHeatmap(
  Matrix = moduleTraitCor,
  
  xLabels = colnames(traitData),
  
  yLabels = names(MEs),
  
  ySymbols = names(MEs),
  
  colorLabels = FALSE,
  
  colors = blueWhiteRed(50),
  
  textMatrix = textMatrix,
  
  setStdMargins = FALSE,
  
  cex.text = 0.8,
  
  cex.lab = 1.0,
  
  cex.axis = 1.0,
  
  zlim = c(-1,1),
  
  main = "Module–Trait Relationships"
)

dev.off()
# =========================
# 26. IDENTIFY BEST MODULE
# =========================

modNames <- substring(names(MEs), 3)

moduleCorrelations <- moduleTraitCor[, "Tumor"]

bestModule <- modNames[
  which.max(abs(moduleCorrelations))
]

bestModule

# =========================
# 27. EXTRACT MODULE GENES
# =========================

moduleGenes <- moduleColors == bestModule

module_gene_names <- colnames(datExpr_filtered)[moduleGenes]

write.csv(
  data.frame(Gene = module_gene_names),
  paste0(bestModule, "_module_genes.csv"),
  row.names = FALSE
)

# =========================
# 28. GENE SIGNIFICANCE
# =========================

GS <- as.data.frame(
  cor(datExpr_filtered,
      traitData$Tumor,
      use = "p")
)

GS_p <- as.data.frame(
  corPvalueStudent(
    as.matrix(GS),
    nrow(datExpr_filtered)
  )
)

# =========================
# 29. MODULE MEMBERSHIP
# =========================

MM <- as.data.frame(
  cor(datExpr_filtered,
      MEs,
      use = "p")
)

MM_p <- as.data.frame(
  corPvalueStudent(
    as.matrix(MM),
    nrow(datExpr_filtered)
  )
)

# =========================
# 30. HUB GENE ANALYSIS
# =========================

hub_df <- data.frame(
  Gene = colnames(datExpr_filtered)[moduleGenes],
  
  GS = GS[moduleGenes, 1],
  
  GS_p = GS_p[moduleGenes, 1],
  
  MM = MM[moduleGenes,
          paste0("ME", bestModule)],
  
  MM_p = MM_p[moduleGenes,
              paste0("ME", bestModule)]
)

# =========================
# 31. FILTER HUB GENES
# =========================

hub_genes <- hub_df[
  abs(hub_df$GS) > 0.3 &
    abs(hub_df$MM) > 0.7 &
    hub_df$GS_p < 0.05 &
    hub_df$MM_p < 0.05,
]

dim(hub_genes) #290   5

head(hub_genes)

# =========================
# 32. SAVE HUB GENES
# =========================

write.csv(
  hub_genes,
  paste0(bestModule, "_HubGenes.csv"),
  row.names = FALSE
)

# =========================
# 33. GS vs MM PLOT
# =========================

MM_selected <- MM[moduleGenes,
                  paste0("ME", bestModule)]

GS_selected <- GS[moduleGenes, 1]

cor_value <- cor(abs(MM_selected),
                 abs(GS_selected))

p_value <- corPvalueStudent(
  cor_value,
  sum(moduleGenes)
)

pdf("6_GS_vs_MM.pdf",
    width = 8,
    height = 6)

plot(
  abs(MM_selected),
  abs(GS_selected),
  xlab = paste("Module Membership in",
               bestModule,
               "module"),
  ylab = "Gene Significance for Tumor",
  main = paste(
    "MM vs GS\ncor =",
    signif(cor_value, 3),
    "\nP =",
    signif(p_value, 3)
  ),
  col = bestModule,
  pch = 16
)

dev.off()

# =========================
# 34. SAVE SESSION
# =========================

save.image("Final_WGCNA_Analysis.RData")

# =========================================================
#                    END OF PIPELINE
# =========================================================