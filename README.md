# 中研院資訊所暑期實習

## Project Overview

- Aim 1:針對深度學習模型已經找到的gene set進行後續分析，比方說enrichment等等。
- Aim 2:透過已有的gene expression data以及某些已知的cell label，預測未知的cell lebel

---

## Current Progress (2026-08-12)

### Completed

- 完成08/12與老師開會的簡報
- 完成RF immune gene的enrichment分析(GSEA only)
- 初步探勘unknown cell的標籤結果
- 分出RF gene list的每個感染病毒子集(WGCNA)

### Current Findings

- 細胞是否感染與細胞類別在卡方檢定下達到顯著(高可信度細胞與預測細胞趨勢不同)
- GCN的表現穩定且非常好
- 細胞存在病毒轉錄體的情況下(不論高低)，高機率被歸類為感染
- 算gene module與感染病毒的關係，需要global sample dataset
- 8個virus family的最佳power不同

---

## Next Steps

- [ ] 想清楚要如何使用WGCNA
- [ ] 將healty的baseline改成1 dpi且vt = 0 的cell 
- [ ] 找出替代scla free 假設的指標
- [ ] 根據RF的統整結果給出生物學意義
- [x] 檢查模型針對unknown cell的標籤


