# 匯入資料
datExpr <- read.table(
  "Alphaviridae_RF_SHAP_expression.tsv",
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)

# WGCNA
library(WGCNA)

options(stringsAsFactors = FALSE)

# 確認資料品質
gsg <- goodSamplesGenes(
  datExpr,
  verbose = 3
)

gsg$allOK


# # 查看sample outlier
# sampleTree <- hclust(
#   dist(datExpr),
#   method = "average"
# )
# 
# plot(
#   sampleTree,
#   main = "Sample clustering to detect outliers",
#   sub = "",
#   xlab = ""
# )

# 轉換expression data
datExpr <- log2(
  datExpr + 1
)

# 選 soft threshold power
powers <- c(
  1:20
)

sft <- pickSoftThreshold(
  datExpr,
  powerVector = powers,
  verbose = 5
)

par(mfrow=c(1,2))

plot(
  sft$fitIndices[,1],
  -sign(sft$fitIndices[,3])*
    sft$fitIndices[,2],
  xlab="Soft Threshold (power)",
  ylab="Scale Free Topology Model Fit",
  type="n"
)

text(
  sft$fitIndices[,1],
  -sign(sft$fitIndices[,3])*
    sft$fitIndices[,2],
  labels=powers
)

abline(
  h=0.8,
  col="red"
)


plot(
  sft$fitIndices[,1],
  sft$fitIndices[,5],
  xlab="Soft Threshold",
  ylab="Mean Connectivity",
  type="n"
)

text(
  sft$fitIndices[,1],
  sft$fitIndices[,5],
  labels=powers
)

# Soft power的選擇
softPower <- 14

# 建立 adjacency matrix
adjacency <- adjacency(
  datExpr,
  power = softPower
)


# 建立 TOM matrix
TOM <- TOMsimilarity(
  adjacency
)

dissTOM <- 1 - TOM


# 建立 gene clustering tree
geneTree <- hclust(
  as.dist(dissTOM),
  method = "average"
)

# # 匯出圖片
# pdf(
#   "geneTree_TOM_clustering.pdf",
#   width = 12,
#   height = 8
# )
# 
# plot(
#   geneTree,
#   main = "Gene clustering on TOM-based dissimilarity",
#   xlab = "",
#   sub = "",
#   labels = FALSE
# )
# 
# dev.off()

dynamicMods <- cutreeDynamic(
  dendro = geneTree,
  distM = dissTOM,
  deepSplit = 2,
  pamRespectsDendro = FALSE,
  minClusterSize = 30
)

dynamicColors <- labels2colors(dynamicMods)

# 計算 Module Eigengenes
MEList <- moduleEigengenes(
  datExpr,
  colors = dynamicColors
)

MEs <- MEList$eigengenes











