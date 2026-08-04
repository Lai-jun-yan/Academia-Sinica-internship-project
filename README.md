# 中研院資訊所暑期實習

## Project Overview

- Aim 1:針對深度學習模型已經找到的gene set進行後續分析，比方說enrichment等等。
- Aim 2:透過已有的gene expression data以及某些已知的cell label，預測未知的cell lebel

---

## Current Progress (2026-08-04)

### Completed

- 使用套件訓練GCN
- 完成RF immune gene的enrichment分析(GSEA only)

### Current Findings

- 細胞是否感染與細胞類別在卡方檢定下達到顯著
- GCN的表現穩定且非常好
- XGBoost所依賴的gene set太小，會造成GSEA結果的解讀有問題

---

## Next Steps

- [ ] GCN與logistic regression和MLP做比較，確定graph是否有幫助
- [ ] 檢查模型針對unknown cell的標籤
- [x] 探討ORA的backgrourd問題
- [ ] 將immune gene資料集跑完XGBoost的enrichment分析
- [ ] 完成第一次針對Aim2的進度報告

