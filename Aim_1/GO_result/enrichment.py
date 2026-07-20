# Pathway database
db = "GO_Biological_Process_2023"

# 先將想分析的欄位取出來，當作enrichment的ranking，之後可以調整
col = "RF_mean_abs_Score"

# 資料前處理

import pandas as pd
data = pd.read_csv(r"C:\Users\USER\Desktop\資訊所實習\計畫\Gene list\consensus_results.Retroviridae.txt",sep="\t")

from sklearn.preprocessing import normalize

# 標準化後算兩種機器學習的Shapley value的平均
data["RF_mean_nor"] = normalize(
    data[["RF_mean_abs_Score"]],
    axis=0
)

data["XGB_mean_nor"] = normalize(
    data[["XGB_mean_abs_Score"]],
    axis=0
)

data["AVG_mean"] = (data["XGB_mean_nor"] + data["RF_mean_nor"])/2

rf_rank = data[["Gene",col]]

# 很多基因的數值為0

rf_rank = rf_rank.loc[rf_rank[col] != 0,:]

rf_rank.columns = ["gene", "score"] # 換col名稱 不確定需不需要

rf_rank = rf_rank.sort_values(
    by="score",
    ascending=False
)

# 取出排名前50的基因
rf_rank_top_100 = rf_rank.iloc[0:50,0:2] 

rf_gene_list = rf_rank_top_100["gene"].tolist() # 可以做OSA

# --------------------------------------
# 跑enrichment
import gseapy as gp

# 先做OSA
res = gp.enrichr(
    gene_list=rf_gene_list, 
    gene_sets=db
)

# 再做GSEA
res_2 = gp.prerank(
    rnk=rf_rank,
    gene_sets=db,
    permutation_num=1000,
    seed=123,
    outdir=None
)

# 結果
print("OSA的結果:")
print(res.results.sort_values(
    by="Adjusted P-value",
    ascending=True
).head(20))
print("-------------------------------------------------------------------------------------------------------------------------")
print("GSEA的結果:")
print(res_2.res2d.sort_values(
    by="FDR q-val",
    ascending=True
).head(20))



