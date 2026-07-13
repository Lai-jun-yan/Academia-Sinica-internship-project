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
# row為 observation (gene)、col為想要分出module的變數
datExpr <- t(expr)

# 使用多核心進行運算
enableWGCNAThreads()

# 檢查資料品質
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

abline(h=0.8,col="red") # 為了讓網絡符合生物學意義，選定R^2當作標準
                        # 但其實如果是cell module，則不一定需要

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

# 算cell之間的關聯性
adjacency <- adjacency(
  datExpr,
  power = softPower,
  type = "signed"
)

# 若靠慮到細胞之間有相連著許多共同的細胞，
# 則會相應的調整細胞之間的相關性
TOM <- TOMsimilarity(
  adjacency,
  TOMType = "signed"
)

# 將相關性轉換成距離，越相關則距離就越近
dissTOM <- 1 - TOM

cellTree <- hclust(
  as.dist(dissTOM),
  method = "average" # cluster之間計算距離的方式
)

# 分出module
dynamicMods <- cutreeDynamic(
  dendro = cellTree,
  distM = dissTOM,
  deepSplit = 4,
  pamRespectsDendro = FALSE,
  minClusterSize = 30
)

dynamicColors <- labels2colors(dynamicMods)

# 不使用 cell type label 的情況下，gene expression 是否包含足夠資訊，
# 讓 network clustering 自動恢復 cell identity？
module_table <- table(
  dynamicColors,
  cell_type
)

# 視覺化不同cell type在各module中的分布情形
# 評估module與cell type之間的對應關係
# 用比例畫圖
pheatmap(
  prop.table(module_table, margin=2), # margin = 1表示module的組成
  cluster_rows = FALSE,               # margin = 2表示cell type的分配
  cluster_cols = TRUE
)

# 將分支圖與module放在一起看
plotDendroAndColors(
  cellTree,
  dynamicColors,
  "Cell Modules",
  dendroLabels = FALSE,
  hang = 0.03
)

###計算ARI，量化結果
module_df <- data.frame(
  cell = colnames(datExpr),
  module = dynamicColors,
  cell_type = cell_type
)

adjustedRandIndex(
  module_df$module,
  module_df$cell_type
)







