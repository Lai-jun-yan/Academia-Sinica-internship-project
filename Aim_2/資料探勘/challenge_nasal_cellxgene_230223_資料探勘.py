import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

adata = sc.read_h5ad(r"C:\Users\USER\Desktop\資訊所實習\計畫\challenge_nasal_cellxgene_230223.h5ad")

# 每個 cell type 保留數量
n_per_type = 100

selected_cells = (
    adata.obs
    .groupby("cell_type", group_keys=False)
    .apply(
        lambda x: x.sample(
            n=min(len(x), n_per_type),
            random_state=123
        )
    )
    .index
)

adata_sub = adata[selected_cells].copy()

adata_hvg = adata_sub[
    :,
    adata_sub.var["vst.variable"]
].copy()

expr_df = pd.DataFrame(
    adata_hvg.X.toarray(),
    index=adata_hvg.obs_names,
    columns=adata_hvg.var_names
)

cell_type_df = adata_hvg.obs[["cell_type"]]

final_df = pd.concat(
    [
        cell_type_df,
        expr_df
    ],
    axis=1
)

print(final_df.head())

final_df["cell_type"].value_counts().plot(
    kind="bar",
    figsize=(10,5)
)

plt.xticks(rotation=90)
plt.show()

final_df = final_df.reset_index()

final_df = final_df.rename(
    columns={"index": "cell_id"}
)

#final_df.to_csv(
#    r"C:\Users\USER\Desktop\資訊所實習\計畫\資料探勘\nasal_scRNA_HVG_100cells.csv",
#    index=False
#)