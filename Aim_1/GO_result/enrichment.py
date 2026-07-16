import pandas as pd

data = pd.read_csv(r"C:\Users\USER\Desktop\資訊所實習\計畫\Gene list\consensus_results.Alphaviridae.txt",sep="\t")

# 先將RF_mean_abs_Score取出來，當作enrichment的ranking

rf_rank = data[["Gene","RF_mean_abs_Score"]]

# 很多基因的數值為0

rf_rank = rf_rank.loc[rf_rank["RF_mean_abs_Score"] != 0,:]

rf_rank.columns = ["gene", "score"]

rf_gene_list = rf_rank["gene"].tolist()

import gseapy as gp

res_2 = gp.prerank(
    rnk=rf_rank,
    gene_sets="KEGG_2021_Human",
    permutation_num=1000,
    seed=123
)

gsea_table = pd.DataFrame(res_2.res2d)

import matplotlib.pyplot as plt


plot_df = gsea_table.head(15).copy()


plt.figure(figsize=(8,6))

plt.scatter(
    plot_df["NES"],
    plot_df["Term"],
    s=plot_df["Tag %"].str.split("/").str[0].astype(int)*20,
    c=plot_df["NOM p-val"],
    cmap="viridis"
)

plt.xlabel("Normalized Enrichment Score (NES)")
plt.ylabel("")
plt.title("KEGG pathway enrichment of RF SHAP-ranked genes")

plt.colorbar(label="Nominal p-value")

plt.tight_layout()
plt.show()

gsea_df = res_2.res2d.copy()

top20 = (
    gsea_df
    .sort_values("NES", ascending=False)
    .head(20)
)

network_df = []

for _, row in top20.iterrows():

    pathway = row["Term"]

    genes = row["Lead_genes"].split(";")

    for gene in genes:

        network_df.append(
            [pathway, gene]
        )


top5 = (
    gsea_df
    .sort_values("NES", ascending=False)
    .head(5)
)

network_df = pd.DataFrame(
    network_df,
    columns=["Pathway","Gene"]
)

network_df = network_df[
    network_df["Pathway"].isin(
        top5["Term"]
    )
]

import networkx as nx


G = nx.Graph()


for _, row in network_df.iterrows():

    G.add_edge(
        row["Pathway"],
        row["Gene"]
    )

plt.figure(figsize=(12,10))


pos = nx.spring_layout(
    G,
    seed=123
)


nx.draw(
    G,
    pos,
    with_labels=True,
    node_size=800,
    font_size=8
)


plt.title(
    "Leading edge gene-pathway network"
)

plt.show()