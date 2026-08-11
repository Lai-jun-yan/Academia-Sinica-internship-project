# 設定工作路徑
setwd("C:/Users/USER/Desktop/資訊所實習/計畫/Gene list/corrected.merged_raw_counts_GRCh38.p13_NCBI.VF_gtr_49_woControl.tsv/WGCNA")

# 匯入資料
datExpr_alpha <- read.table(
  "Alphaviridae_RF_SHAP_expression.tsv",
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)

datExpr_Coron <- read.table(
  "Coronaviridae_RF_SHAP_expression.tsv",
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)

datExpr_Flavi <- read.table(
  "Flaviviridae_RF_SHAP_expression.tsv",
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)

datExpr_Orthoher <- read.table(
  "Orthoherpesviridae_RF_SHAP_expression.tsv",
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)

datExpr_Orthomy <- read.table(
  "Orthomyxoviridae_RF_SHAP_expression.tsv",
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)

datExpr_Picor <- read.table(
  "Picornaviridae_RF_SHAP_expression.tsv",
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)

datExpr_Reovi <- read.table(
  "Reoviridae_RF_SHAP_expression.tsv",
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)

datExpr_Retro <- read.table(
  "Retroviridae_RF_SHAP_expression.tsv",
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)

# 將8個virus family合併成list
datasets <- list(
  Alpha = datExpr_alpha,
  Coron = datExpr_Coron,
  Flavi = datExpr_Flavi,
  Orthoher = datExpr_Orthoher,
  Orthomy = datExpr_Orthomy,
  Picor = datExpr_Picor,
  Reovi = datExpr_Reovi,
  Retro = datExpr_Retro
)

# WGCNA
library(WGCNA)

# 對每個virus family做QC
gsg_list <- lapply(
  datasets,
  function(x) {
    goodSamplesGenes(
      x,
      verbose = 3
    )
  }
)

# 找出共同品質有過關的基因
good_gene_names <- lapply(
  seq_along(datasets),
  function(i) {
    colnames(datasets[[i]])[gsg_list[[i]]$goodGenes]
  }
)

names(good_gene_names) <- names(datasets)

common_genes <- Reduce(
  intersect,
  good_gene_names
)

length(common_genes)

# 然後統一基因
common_genes <- sort(common_genes)

datasets_common <- lapply(
  datasets,
  function(x) {
    x[, common_genes, drop = FALSE]
  }
)

# log transformation
datasets_common <- lapply(
  datasets_common,
  function(x) {
    log2(x + 1)
  }
)

# 選 soft power
powers <- 1:20

sft_list <- lapply(
  datasets_common,
  function(x) {
    pickSoftThreshold(
      x,
      powerVector = powers,
      verbose = 0
    )
  }
)

sft_table <- do.call( # 整理成表格
  rbind,
  lapply(
    names(sft_list),
    function(v) {
      data.frame(
        Virus = v,
        sft_list[[v]]$fitIndices
      )
    }
  )
)

# 把每個virus最佳的power選出來
softPower_table <- do.call(
  rbind,
  lapply(
    names(sft_list),
    function(v) {
      
      df <- sft_list[[v]]$fitIndices
      
      idx <- which(df$SFT.R.sq >= 0.8)[1]
      
      data.frame(
        Virus = v,
        SelectedPower = ifelse(
          length(idx) == 0,
          NA,
          df$Power[idx]
        )
      )
    }
  )
)

softPower_table

# 畫圖比較

# 存成pdf檔
pdf(
  "WGCNA_soft_threshold.pdf",
  width = 12,
  height = 6
)

# ============================================
# 設定 virus family 顏色
# ============================================

virus_colors <- c(
  Alpha    = "#E41A1C",
  Coron    = "#377EB8",
  Flavi    = "#4DAF4A",
  Orthoher = "#984EA3",
  Orthomy  = "#FF7F00",
  Picor    = "#A65628",
  Reovi    = "#F781BF",
  Retro    = "#17BECF"
)

par(mfrow = c(1, 2))


# ============================================
# 左圖：SFT.R²
# ============================================

plot(
  NA,
  xlim = range(powers),
  ylim = c(0, 1),
  xlab = "Soft Threshold Power",
  ylab = "Scale-Free Topology Fit (SFT.R²)",
  main = "Scale-Free Topology Fit",
  cex.lab = 1.1,
  cex.main = 1.2
)

for (v in names(sft_list)) {
  
  fit <- sft_list[[v]]$fitIndices
  
  lines(
    fit[, "Power"],
    fit[, "SFT.R.sq"],
    type = "b",
    pch = 16,
    col = virus_colors[v],
    lwd = 2
  )
  
}

# R² = 0.8 threshold
abline(
  h = 0.8,
  lty = 2,
  lwd = 1.5
)


# ============================================
# 右圖：Mean Connectivity
# ============================================

all_mean_k <- unlist(
  lapply(
    sft_list,
    function(x) x$fitIndices[, "mean.k."]
  )
)

plot(
  NA,
  xlim = range(powers),
  ylim = range(all_mean_k),
  xlab = "Soft Threshold Power",
  ylab = "Mean Connectivity",
  main = "Mean Connectivity",
  cex.lab = 1.1,
  cex.main = 1.2
)

for (v in names(sft_list)) {
  
  fit <- sft_list[[v]]$fitIndices
  
  lines(
    fit[, "Power"],
    fit[, "mean.k."],
    type = "b",
    pch = 16,
    col = virus_colors[v],
    lwd = 2
  )
  
}


# ============================================
# Legend
# ============================================

legend(
  "topright",
  legend = names(sft_list),
  col = virus_colors[names(sft_list)],
  pch = 16,
  lty = 1,
  lwd = 2,
  cex = 0.75,
  bty = "n"
)


par(mfrow = c(1, 1))

# 關閉圖片裝置，正式儲存
dev.off()



























