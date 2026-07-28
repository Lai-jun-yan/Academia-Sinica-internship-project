import scanpy as sc

adata = sc.read_h5ad(r"C:\Users\USER\Desktop\資訊所實習\計畫\資料探勘\Ravindra2021.raw_count.stdprep.h5ad")

# 整體品質管控
adata_qc = adata[
    (adata.obs["pmito"] < 0.2) # 粒線體的的RNA比例超過0.2，代表細胞品質不好
    &
    (adata.obs["n_genes"] > 500) # 基因表達數量需要高於500(根據分布，已屬於寬鬆的結果)，較具有特徵
].copy()

# 看一下若有偵測到病毒轉錄體的情況下，病毒轉錄體的分布，以決定感染標籤的標準
print("偵測到病毒轉錄體的情況下，病毒轉錄體的分布:")
print(adata_qc.obs[adata_qc.obs["Viral_transcript"] > 0]["Viral_transcript"].describe())
print("---------------------------------------------------")

# 偵測到病毒轉錄體的情況下，病毒轉錄體的分布:
# count    38169.000000
# mean        87.941602
# std        598.887808
# min          1.000000
# 25%          1.000000
# 50%          3.000000
# 75%          4.000000
# max      18682.000000
# 極度右偏，感染標準用百分位數而非平均+標準差

# 由資料分布決定 virus threshold
virus_threshold = (
    adata_qc.obs.loc[
        adata_qc.obs["Viral_transcript"] > 0, # 未偵測到病毒轉錄體的細胞不考慮
        "Viral_transcript"
    ]
    .quantile(0.95) # 取95百分位數當作標準
)

positive_cell = adata_qc[
    (adata_qc.obs["Viral_transcript"] >= virus_threshold)
    &
    (adata_qc.obs["Condition"] != "Mock") # 不能是確定未感染的細胞
    &
    (adata_qc.obs["n_genes"] > 1000) # 高可信標籤需要更嚴格的標準(避免Empty droplet)
]

# 選取未感染的細胞
negative_cell = adata_qc[
    (adata_qc.obs["Viral_transcript"] == 0) # 沒有偵測到病毒轉錄體
    &
    (adata_qc.obs["Condition"] == "Mock") # 確定未感染的細胞
]

# 確認高可信度細胞的資料品質
print("高可信度細胞的資料品質:")
print("positive")
print(positive_cell.obs[
    ["Viral_transcript", "n_genes", "pmito"]
].describe())
print("")
print("negative")
print(negative_cell.obs[
    ["Viral_transcript", "n_genes", "pmito"]
].describe())
print("---------------------------------------------------")

# 高可信度細胞的資料品質:
# positive
#        Viral_transcript      n_genes        pmito
# count       1713.000000  1713.000000  1713.000000
# mean        1694.311734  3824.687682     0.070038
# std         2216.431011  1769.735079     0.031631
# min          165.000000  1009.000000     0.001452
# 25%          425.000000  2327.000000     0.049580
# 50%          872.000000  3649.000000     0.069454
# 75%         2004.000000  5103.000000     0.088599
# max        18682.000000  9009.000000     0.196126

# negative
#        Viral_transcript       n_genes         pmito
# count           22080.0  22080.000000  22080.000000
# mean                0.0   2656.235779      0.043987
# std                 0.0   1353.388096      0.023677
# min                 0.0    501.000000      0.000116
# 25%                 0.0   1643.000000      0.027254
# 50%                 0.0   2272.000000      0.040484
# 75%                 0.0   3383.250000      0.056685
# max                 0.0   8773.000000      0.199321

# 取得已標籤 cell 的 index
positive_idx = positive_cell.obs_names
negative_idx = negative_cell.obs_names

# unknown = QC後資料扣掉positive和negative
unknown_cell = adata_qc[
    ~adata_qc.obs_names.isin(positive_idx.union(negative_idx))
]

# 檢查不確定細胞的資料品質
print("未知狀態細胞的資料品質:")
print(unknown_cell.obs[
    ["Viral_transcript", "n_genes", "pmito"]
].describe())
print("---------------------------------------------------")

#        Viral_transcript       n_genes         pmito
# count      50295.000000  50295.000000  50295.000000
# mean           9.032449   3120.482593      0.080332
# std          113.933869   1671.276747      0.030199
# min            0.000000    501.000000      0.000286
# 25%            0.000000   1813.000000      0.060610
# 50%            2.000000   2842.000000      0.078488
# 75%            3.000000   4066.000000      0.097949
# max         8650.000000  10571.000000      0.199936

# 在資料中新增標籤欄位
# 初始化 unknown
adata_qc.obs["label"] = -1


# infected
adata_qc.obs.loc[
    positive_cell.obs_names,
    "label"
] = 1


# uninfected
adata_qc.obs.loc[
    negative_cell.obs_names,
    "label"
] = 0

adata_qc.obs["label_source"] = "unknown"

adata_qc.obs.loc[
    positive_cell.obs_names,
    "label_source"
] = "high_confidence_infected"

adata_qc.obs.loc[
    negative_cell.obs_names,
    "label_source"
] = "high_confidence_uninfected"


import pandas as pd

print("看一下細胞類型與感染狀態是否有關:")
print(pd.crosstab(
    adata_qc.obs["ctype"],
    adata_qc.obs["label"]
))

# 看一下細胞類型與感染狀態是否有關:
# label                    -1     0     1
# ctype                                  
# BC/Club                2351  1521   124
# Basal cells           19330  9549   169
# Ciliated cells        15172  6208  1222
# Club cells            12302  4493   194
# Goblet cells             74    43     0
# Ionocytes               530   123     2
# Neuroendocrine cells    299    74     2
# Tuft cells              237    69     0
# 後續評估GCN模型表現，可以根據cell type stratification






