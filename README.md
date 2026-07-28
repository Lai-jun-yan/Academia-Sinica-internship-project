# 中研院資訊所暑期實習

## Project Overview

- Aim 1:針對深度學習模型已經找到的gene set進行後續分析，比方說enrichment等等。
- Aim 2:透過已有的gene expression data以及某些已知的cell label，預測未知的cell lebel

---

## Current Progress (2026-07-28)

### Completed

- 完成初步的Ravindra2021.raw_count.stdprep.h5ad資料探勘
- 建立選擇高可信度細胞標籤的框架

### Current Findings

- 病毒轉錄體的分布呈現右偏態
- 細胞是否感染與細胞類別可能有關

---

## Next Steps

- [ ] 針對Ravindra2021.raw_count.stdprep.h5ad的基因做特徵選取
- [ ] 探討ORA的backgrourd問題
- [ ] 將目前取得的gene list跑完enrcihment，用Mean的score應該就可以了
- [x] 透過 Ravindra2021.raw_count.stdprep.h5ad 分出判別細胞感染的threshold

