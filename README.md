# 中研院資訊所暑期實習

## Project Overview

- Aim 1:針對深度學習模型已經找到的gene set進行後續分析，比方說enrichment等等。
- Aim 2:透過已有的gene expression data以及某些已知的cell label，預測未知的cell lebel

---

## Current Progress (2026-08-07)

### Completed

- 使用套件訓練GCN
- 完成RF immune gene的enrichment分析(GSEA only)
- 初步探勘unknown cell的標籤結果
- 分出RF gene list的每個感染病毒子集(WGCNA)

### Current Findings

- 細胞是否感染與細胞類別在卡方檢定下達到顯著
- GCN的表現穩定且非常好
- XGBoost所依賴的gene set太小，會造成GSEA結果的解讀有問題
- 細胞存在病毒轉錄體的情況下(不論高低)，高機率被歸類為感染
- 算gene module與感染病毒的關係，需要global sample dataset

---

## Next Steps

- [ ] 每個病毒子集跑WGCNA
- [ ] 合併每個病毒子集(global)
- [x] 檢查模型針對unknown cell的標籤
- [x] 探討ORA的backgrourd問題
- [ ] 探討RF可能學習到與哪些生物意義相關


