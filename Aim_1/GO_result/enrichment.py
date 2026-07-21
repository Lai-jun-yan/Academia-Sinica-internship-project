# Pathway database
db = "WikiPathways_2024_Human"

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

# OSA顯著的結果
ora_df = res.results.copy()


# 取顯著 pathway
sig_df = ora_df[
    ora_df["Adjusted P-value"] < 0.05
].copy()


# 依顯著性排序
sig_df = sig_df.sort_values(
    "Adjusted P-value",
    ascending=True
)

# GSEA顯著的結果
sig_gsea = res_2.res2d[
    res_2.res2d["FDR q-val"] < 0.25
]

sig_gsea = sig_gsea.sort_values(
    by="FDR q-val",
    ascending=True
)

# ====================================================
# ORA與GSEA共同顯著的pathway
# ====================================================

# ORA顯著pathway
ora_sig = sig_df.copy()

# GSEA顯著pathway
gsea_sig = sig_gsea.copy()

# 找共同pathway
common_pathways = sorted(
    set(ora_sig["Term"]).intersection(
        set(gsea_sig["Term"])
    )
)

print("共同顯著Pathway數量:", len(common_pathways))

# ====================================================
# 建立整理表
# ====================================================

if len(common_pathways) == 0:

    print("沒有ORA與GSEA同時顯著的Pathway")
    
    common_table = pd.DataFrame(
        columns=[
            "Pathway",
            "ORA_FDR",
            "GSEA_FDR",
            "ORA_Genes",
            "Leading_edge",
            "Union_genes",
            "Union_gene_count"
        ]
    )

else:

    common_table = []

    for pathway in common_pathways:

        ora_row = ora_sig.loc[
            ora_sig["Term"] == pathway
        ].iloc[0]

        gsea_row = gsea_sig.loc[
            gsea_sig["Term"] == pathway
        ].iloc[0]


        ora_genes = set(
            ora_row["Genes"].split(";")
        )

        lead_genes = set(
            gsea_row["Lead_genes"].split(";")
        )


        union_genes = sorted(
            ora_genes.union(lead_genes)
        )


        common_table.append({
            "Pathway": pathway,
            "ORA_FDR": ora_row["Adjusted P-value"],
            "GSEA_FDR": gsea_row["FDR q-val"],
            "ORA_Genes": ";".join(sorted(ora_genes)),
            "Leading_edge": ";".join(sorted(lead_genes)),
            "Union_genes": ";".join(union_genes),
            "Union_gene_count": len(union_genes)
        })


    common_table = pd.DataFrame(common_table)

    common_table = common_table.sort_values(
        by=["ORA_FDR","GSEA_FDR"]
    )


print(common_table)