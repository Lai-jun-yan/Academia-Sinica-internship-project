# 中研院資訊所暑期實習

## Project Overview

- Aim 1:針對深度學習模型已經找到的gene set進行後續分析，比方說enrichment等等。
- Aim 2:透過已有的gene expression data以及某些已知的cell label，預測未知的cell lebel

---

## Current Progress (2026-07-23)

### Completed

- 完成RF mean的database分析(KEGG、GO、Reactome、MSigDB_Hallmark、WikiPathways)
- 下載Ravindra2021.raw_count.stdprep.h5ad的資料
- 大致看完 Human SARS-CoV-2 challenge uncovers local and systemic response dynamics

### Current Findings

- RF相關的結果很多，需要有系統的分析
- 針對RF mean的部分，只有Retroviridae跟Picornaviridae在OSA以及GSEA同時擁有共同顯著的pathway
- +-標準差去找出高可信度的細胞標籤

---

## Next Steps

- [x] 與Serina討論想法
- [x] 解決GSEA gene size太小的問題
- [x] 針對RF mean的結果，初步給出生物學的故事
- [ ] 將目前取得的gene list跑完enrcihment，用Mean的score應該就可以了
- [ ] 透過 Ravindra2021.raw_count.stdprep.h5ad 分出判別細胞感染的threshold

