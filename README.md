# 中研院資訊所暑期實習

## Project Overview

- Aim 1:針對深度學習模型已經找到的gene set進行後續分析，比方說enrichment等等。
- Aim 2:透過已有的gene expression data以及某些已知的cell label，預測未知的cell lebel

---

## Current Progress (2026-07-14)

### Completed

- 大致看過暑期實習計畫的參考文獻
- 看了enrichment的review paper，熟悉概念
- 將challenge_nasal_cellxgene_230223.h5ad的cell expression profile部分資料拿出來跑WGCNA
- 用challenge_nasal_cellxgene_230223.h5ad取出部分細胞做KNN，接著嘗試跑GCN

### Current Findings

- 對於WGCNA的應用，我覺得兩個實習目標都適用，但目標二需要測試效果
- 如果用WGCNA找cell module，目前測試的效果不太好，似乎不好應用在目標二方面
- 如果用KNN的方式畫出cell graph，再用GCN，預測效果看起來不錯

---

## Next Steps

- [x] 完成文獻回顧的簡報
- [ ] 與Serina討論想法
- [x] 確認實際的資料型態
- [ ] 確認KNN的過程正確
- [ ] 取得gene list，開始針對目標一進行enrichment
- [x] 將試跑WGCNA的過程與簡報統整出來
