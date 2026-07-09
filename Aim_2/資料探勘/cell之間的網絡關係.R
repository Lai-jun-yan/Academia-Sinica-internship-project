library(WGCNA)
library(mclust)
library(pheatmap)

# 讀檔
expr <- read.csv(
  "nasal_scRNA_HVG_100cells.csv",
  row.names = 1,
  check.names = FALSE
)

# 將cell type欄位先分出來
cell_type <- expr$cell_type

expr <- expr[, colnames(expr) != "cell_type"]

# 根據WGCNA的資料格式進行轉至
datExpr <- t(expr)

# 檢查資料品質
enableWGCNAThreads()

gsg <- goodSamplesGenes(
  datExpr,
  verbose = 3
)

gsg$allOK


# 選定相關性的指數
powers <- c(1:10, seq(12,20,2))

sft <- pickSoftThreshold(
  datExpr,
  powerVector = powers,
  networkType = "signed",
  verbose = 5
)

par(mfrow=c(1,2))

plot(
  sft$fitIndices[,1],
  -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
  xlab="Soft Threshold",
  ylab="Scale Free Topology Model Fit",
  type="n"
)

text(
  sft$fitIndices[,1],
  -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
  labels=powers,
  col="red"
)

abline(h=0.8,col="red")

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
  labels=powers,
  col="red"
)


# 建立 adjacency matrix
softPower <- 14

adjacency <- adjacency(
  datExpr,
  power = softPower,
  type = "signed"
)

TOM <- TOMsimilarity(
  adjacency,
  TOMType = "signed"
)

dissTOM <- 1 - TOM

cellTree <- hclust(
  as.dist(dissTOM),
  method = "average"
)

dynamicMods <- cutreeDynamic(
  dendro = cellTree,
  distM = dissTOM,
  deepSplit = 4,
  pamRespectsDendro = FALSE,
  minClusterSize = 30
)

dynamicColors <- labels2colors(dynamicMods)

plotDendroAndColors(
  cellTree,
  dynamicColors,
  "Cell Modules",
  dendroLabels = FALSE,
  hang = 0.03
)

# 判斷module與cell type是否一致
module_table <- table(
  dynamicColors,
  cell_type
)

module_table

###
module_df <- data.frame(
  cell = colnames(datExpr),
  module = dynamicColors,
  cell_type = cell_type
)

head(module_df)

module_prop <- prop.table(
  table(
    module_df$module,
    module_df$cell_type
  ),
  margin = 1
)

module_prop

pheatmap(
  module_prop,
  cluster_rows = FALSE,
  cluster_cols = TRUE
)

adjustedRandIndex(
  module_df$module,
  module_df$cell_type
)







