# 中研院資訊所暑期實習

## Project Overview

- Aim 1:針對深度學習模型已經找到的gene set進行後續分析，比方說enrichment等等。
- Aim 2:透過已有的gene expression data以及某些已知的cell label，預測未知的cell lebel

---

## Current Progress (2026-08-17)

### Completed

- 完成08/12與老師開會的簡報
- 完成RF immune gene的enrichment分析(GSEA only)
- Aim 2跑logistic regression
- 以1dpi當作baseline，重新跑GCN
- 初步探討GCN與Logistic的差別

### Current Findings

- 細胞是否感染與細胞類別在卡方檢定下達到顯著(高可信度細胞與預測細胞趨勢不同)
- Logistic regression的表現效果比GCN好
- 將mock剔除後，GCN的表現下降，但還是不錯
- GCN vs Logistic對於unknown cell的預測有差別
- GCN & Logistic仍然有病毒轉錄體高，而感染機率上升的趨勢

---

## Next Steps

- [ ] 想清楚要如何使用WGCNA
- [x] 將healty的baseline改成1 dpi且vt = 0 的cell 
- [ ] 找出替代scla free 假設的指標
- [ ] 根據RF的統整結果給出生物學意義
- [x] 使用統計模型(如羅吉斯回歸)與GCN比較


