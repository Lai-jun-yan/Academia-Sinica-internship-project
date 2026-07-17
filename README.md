# 中研院資訊所暑期實習

## Project Overview

- Aim 1:針對深度學習模型已經找到的gene set進行後續分析，比方說enrichment等等。
- Aim 2:透過已有的gene expression data以及某些已知的cell label，預測未知的cell lebel

---

## Current Progress (2026-07-17)

### Completed

- 看了enrichment的review paper，熟悉概念
- 完成OSA跟GSEA的pipeline
- 將challenge_nasal_cellxgene_230223.h5ad的cell expression profile部分資料拿出來跑WGCNA
- 用challenge_nasal_cellxgene_230223.h5ad取出部分細胞做KNN，接著嘗試跑GCN

### Current Findings

- GSEA的input gene數量不能太少，不然會無法分析
- 如果用WGCNA找cell module，目前測試的效果不太好，似乎不好應用在目標二方面
- 如果用KNN的方式畫出cell graph，再用GCN，預測效果看起來不錯

---

## Next Steps

- [x] 與Serina討論想法
- [ ] 解決GSEA gene size太小的問題
- [ ] 將目前取得的gene list跑完enrcihment，用AVG的score應該就可以了

