# 中研院資訊所暑期實習

## Project Overview

- Aim 1:針對深度學習模型已經找到的gene set進行後續分析，比方說enrichment等等。
- Aim 2:透過已有的gene expression data以及某些已知的cell label，預測未知的cell lebel

---

## Current Progress (2026-07-30)

### Completed

- 建立選擇高可信度細胞標籤的框架
- 匯出訓練GCN所需要的資料
- 使用套件訓練GCN

### Current Findings

- 病毒轉錄體的分布呈現右偏態
- 細胞是否感染與細胞類別在卡方檢定下達到顯著
- GCN的表現穩定且非常好

---

## Next Steps

- [ ] GCN與logistic regression和MLP做比較，確定graph是否有幫助
- [ ] 檢查模型針對unknown cell的標籤
- [ ] 探討ORA的backgrourd問題
- [ ] 將目前取得的gene list跑完enrcihment，用Mean的score應該就可以了
- [ ] 完成第一次針對Aim2的進度報告

