# Pathway database
db = "KEGG_2021_Human"

# 先將想分析的欄位取出來，當作enrichment的ranking，之後可以調整
col = "RF_mean_abs_Score"


# 資料前處理

import pandas as pd
data = pd.read_csv(r"C:\Users\USER\Desktop\資訊所實習\計畫\Gene list\consensus_results.Alphaviridae.txt",sep="\t")

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

# 取出排名前100的基因
rf_rank_top_100 = rf_rank.iloc[0:100,0:2] 

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

# ----------------------------------------

# 對機器學習模型最重要的前20個Gene
top20_gene = (
    rf_rank.head(20).copy()
)

import matplotlib.pyplot as plt

plt.figure(figsize=(8,6))

plt.barh(
    top20_gene["gene"][::-1],
    top20_gene["score"][::-1]
)

plt.xlabel("RF SHAP importance")
plt.ylabel("Gene")

plt.title(
    "Top 20 genes ranked by RF SHAP importance"
)

plt.tight_layout()

plt.show()

# Shapley value distribution
shap_rank = (
    data["RF_mean_abs_Score"]
    .sort_values(
        ascending=False
    )
    .reset_index(drop=True)
)


plt.figure(figsize=(8,5))


plt.plot(
    shap_rank
)


plt.xlabel(
    "Gene rank"
)

plt.ylabel(
    "RF SHAP importance"
)


plt.title(
    "Ranked distribution of RF SHAP importance"
)


plt.tight_layout()

plt.show()

# OSA前10個顯著的pathway
ora = res.results.copy()

ora[["Gene_count","Pathway_size"]] = (
    ora["Overlap"]
    .str.split("/", expand=True)
    .astype(int)
)


ora["Gene_ratio"] = (
    ora["Gene_count"] /
    ora["Pathway_size"]
)

ora_top10 = (
    ora
    .sort_values(
        "Adjusted P-value"
    )
    .head(10)
)

plt.figure(figsize=(8,6))


plt.scatter(
    ora_top10["Gene_ratio"],
    ora_top10["Term"],
    s=ora_top10["Gene_count"]*30, # leading edge的基因數量越大，圓越大
    c=ora_top10["Adjusted P-value"],
    cmap="viridis"
)


plt.xlabel("Gene ratio")

plt.ylabel("Pathway")

plt.title(
    "Top enriched pathways among top SHAP genes"
)


plt.colorbar(
    label="Adjusted P-value"
)


plt.tight_layout()

plt.show()

# GSEA依照NES前10名以及後10名排列
gsea = res_2.res2d.copy()

gsea_top = pd.concat([
    gsea.sort_values("NES", ascending=False).head(10),
    gsea.sort_values("NES", ascending=True).head(10)
])

gsea_top["Gene_count"] = (
    gsea_top["Tag %"]
    .str.split("/")
    .str[0]
    .astype(int)
)

plt.figure(figsize=(8,6))


plt.scatter(
    gsea_top["NES"],
    gsea_top["Term"],
    s=gsea_top["Gene_count"]*20,
    c=gsea_top["FWER p-val"],
    cmap="viridis"
)

plt.axvline(
    0,
    linestyle="--"
)

plt.xlim(
    gsea_top["NES"].min()-0.1,
    gsea_top["NES"].max()+0.1
)

plt.xlabel(
    "Normalized Enrichment Score (NES)"
)

plt.ylabel(
    "Pathway"
)


plt.title(
    "Exploratory GSEA of SHAP-ranked genes"
)


plt.colorbar(
    label="FWER p-val"
)


plt.tight_layout()

plt.show()

# OSA跟GSEA放在一起做比較
ora_compare = ora[
    [
        "Term",
        "Adjusted P-value"
    ]
].copy()


ora_compare.columns = [
    "Pathway",
    "OSA_FDR"
]

gsea_compare = gsea[
    [
        "Term",
        "NES",
        "FDR q-val",
        "Tag %"
    ]
].copy()


gsea_compare.columns = [
    "Pathway",
    "GSEA_NES",
    "GSEA_FDR",
    "Leading_edge"
]

compare_table = pd.merge(
    ora_compare,
    gsea_compare,
    on="Pathway",
    how="inner"
)

compare_table[["Leading_gene","Pathway_size"]] = (
    compare_table["Leading_edge"]
    .str.split("/", expand=True)
    .astype(int)
)

compare_table["Leading_ratio"] = (
    compare_table["Leading_gene"] /
    compare_table["Pathway_size"]
)

compare_table.sort_values(
    "GSEA_NES",
    ascending=False
).head(20)

import numpy as np

heat = (
    compare_table
    .sort_values(
        "GSEA_NES",
        ascending=False
    )
    .copy()
)

heat["ORA"] = -np.log10(
    heat["OSA_FDR"]
)

heat["GSEA"] = heat["GSEA_NES"]

heat["Leading"] = heat["Leading_ratio"]

heat = heat[
    [
        "Pathway",
        "ORA",
        "GSEA",
        "Leading"
    ]
]

heat = heat.set_index(
    "Pathway"
)

from sklearn.preprocessing import StandardScaler

scaler = StandardScaler()

heat_scaled = pd.DataFrame(
    scaler.fit_transform(heat),
    index=heat.index,
    columns=heat.columns
)

plt.figure(figsize=(6,10))

plt.imshow(
    heat_scaled,
    aspect="auto",
    cmap="RdBu_r" # 轉成-log10
)

# x 軸
plt.xticks(
    range(len(heat.columns)),
    heat.columns
)

# y 軸
plt.yticks(
    range(len(heat.index)),
    heat.index,
    fontsize=8
)

plt.colorbar(label="Value")

plt.title("Comparison of ORA and GSEA Results")

plt.tight_layout()

plt.show()

# 畫關係網絡圖
import networkx as nx
import matplotlib.pyplot as plt


##########################################
# 取 GSEA NES 前10個共同 pathway
##########################################

top_pathway = (
    compare_table # 兩種結果共有的pathway
    .sort_values(
        "GSEA_NES",
        ascending=False
    )
    .head(10)["Pathway"]
    .tolist()
)


network_df = gsea[
    gsea["Term"].isin(top_pathway)
].copy()


##########################################
# Gene SHAP importance
##########################################

shap_dict = dict(
    zip(
        data["Gene"],
        data["RF_mean_abs_Score"]
    )
)


##########################################
# Pathway NES
##########################################

nes_dict = dict(
    zip(
        network_df["Term"],
        network_df["NES"]
    )
)


##########################################
# 建立 network
##########################################

G = nx.Graph()


for _, row in network_df.iterrows():

    pathway = row["Term"]

    genes = row["Lead_genes"].split(";")

    for gene in genes:

        G.add_edge(
            gene,
            pathway
        )

##########################################
# 列出使用到的 pathway
##########################################

pathway_nodes = top_pathway

print("Pathways used in network:\n")

for pathway in pathway_nodes:
    print(pathway)


##########################################
# 分類 node
##########################################

gene_nodes = [
    node for node in G.nodes()
    if node in shap_dict
]


pathway_nodes = [
    node for node in G.nodes()
    if node not in shap_dict
]

##########################################
# 計算 gene degree
##########################################

gene_degree = []

for gene in gene_nodes:

    gene_degree.append(
        [
            gene,
            G.degree(gene),
            shap_dict[gene]
        ]
    )


hub_table = pd.DataFrame(
    gene_degree,
    columns=[
        "Gene",
        "Pathway_count",
        "SHAP_importance"
    ]
)


hub_table = (
    hub_table
    .sort_values(
        by=[
            "Pathway_count",
            "SHAP_importance"
        ],
        ascending=False
    )
)


print(
    hub_table.head(10)
)

##########################################
# node size
##########################################

gene_sizes = [
    shap_dict[gene]*20000 + 100
    for gene in gene_nodes
]


pathway_sizes = [
    abs(nes_dict[pathway])*3000
    for pathway in pathway_nodes
]


##########################################
# 畫圖
##########################################

plt.figure(figsize=(12,10))


pos = nx.spring_layout(
    G,
    seed=123,
    k=1
)


# pathway node

nx.draw_networkx_nodes(
    G,
    pos,
    nodelist=pathway_nodes,
    node_color="salmon",
    node_size=pathway_sizes
)


# gene node

nx.draw_networkx_nodes(
    G,
    pos,
    nodelist=gene_nodes,
    node_color="skyblue",
    node_size=gene_sizes
)


# edges

nx.draw_networkx_edges(
    G,
    pos,
    alpha=0.5
)


# labels

nx.draw_networkx_labels(
    G,
    pos,
    font_size=9
)


plt.title(
    "Top 10 GSEA Pathways and Leading-edge Genes"
)


plt.axis("off")

plt.tight_layout()

plt.show()