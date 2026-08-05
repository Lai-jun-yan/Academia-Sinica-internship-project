import numpy as np
import pandas as pd
from scipy import sparse

import torch
from torch_geometric.data import Data


# ======================
# Load data
# ======================

# node features (cell × gene)
X = np.load(
    r"C:\Users\USER\Desktop\資訊所實習\計畫\資料探勘\Ravindra2021.raw_count.stdprep.h5ad的細胞標籤\GCN_features.npy"
)

# print("Feature shape:", X.shape)


# adjacency matrix
A = sparse.load_npz(
    r"C:\Users\USER\Desktop\資訊所實習\計畫\資料探勘\Ravindra2021.raw_count.stdprep.h5ad的細胞標籤\GCN_adjacency.npz"
)

# metadata
metadata = pd.read_csv(
    r"C:\Users\USER\Desktop\資訊所實習\計畫\資料探勘\Ravindra2021.raw_count.stdprep.h5ad的細胞標籤\GCN_metadata.csv",
    index_col=0
)

# print(X.shape)
# print(A.shape)
# print(metadata.shape)

labels = metadata["label"].values

y = torch.tensor(
    labels,
    dtype=torch.long
)

from sklearn.model_selection import train_test_split


# 有標籤的位置
label_idx = np.where(labels != -1)[0]


# train / validation
train_idx, val_idx = train_test_split(
    label_idx,
    test_size=0.2,
    stratify=labels[label_idx],
    random_state=42
)


train_mask = torch.zeros(
    len(labels),
    dtype=torch.bool
)

val_mask = torch.zeros(
    len(labels),
    dtype=torch.bool
)


train_mask[train_idx] = True
val_mask[val_idx] = True


# print("Train:", train_mask.sum())
# print("Validation:", val_mask.sum())
# print("Unknown:", (labels==-1).sum())


edge_index = torch.tensor(
    np.vstack(A.nonzero()),
    dtype=torch.long
)


# print(edge_index.shape)
data = Data(
    x=torch.tensor(
        X,
        dtype=torch.float
    ),

    edge_index=edge_index,

    y=y,

    train_mask=train_mask,

    val_mask=val_mask
)


# print(data)

import torch.nn.functional as F
from torch_geometric.nn import GCNConv


class GCN(torch.nn.Module):

    def __init__(self):
        super().__init__()

        self.conv1 = GCNConv(
            3000,
            64
        )

        self.conv2 = GCNConv(
            64,
            2
        )


    def forward(self, x, edge_index):

        x = self.conv1(
            x,
            edge_index
        )

        x = F.relu(x)

        x = F.dropout(
            x,
            p=0.5,
            training=self.training
        )


        x = self.conv2(
            x,
            edge_index
        )

        return x

device = torch.device(
    "cuda" if torch.cuda.is_available()
    else "cpu"
)


data = data.to(device)

model = GCN().to(device)


optimizer = torch.optim.Adam(
    model.parameters(),
    lr=0.01,
    weight_decay=5e-4
)


loss_fn = torch.nn.CrossEntropyLoss()

print("以下為訓練過程:")

for epoch in range(100):

    model.train()

    optimizer.zero_grad()


    out = model(
        data.x,
        data.edge_index
    )


    loss = loss_fn(
        out[data.train_mask],
        data.y[data.train_mask]
    )


    loss.backward()

    optimizer.step()



    # validation

    model.eval()

    with torch.no_grad():

        val_out = model(
            data.x,
            data.edge_index
        )

        pred = val_out.argmax(dim=1) # 選擇模型輸出最高的類別

        correct = (
            pred[data.val_mask]
            ==
            data.y[data.val_mask]
        ).sum()

        acc = correct / data.val_mask.sum()


    if epoch % 10 == 0:

        print(
            epoch,
            "Loss:",
            loss.item(),
            "Val Acc:",
            acc.item()
        )

print("----------------------------------------------------------------")

print("模型的表現:")
from sklearn.metrics import confusion_matrix, classification_report


model.eval()

with torch.no_grad():

    out = model(
        data.x,
        data.edge_index
    )

    pred = out.argmax(dim=1) # 選擇模型輸出最高的類別

        # 感染機率 (class 1)
    prob = torch.softmax(
        out,
        dim=1
    )[:,1]


metadata["predict_label"] = pred
metadata["infect_prob"] = prob.cpu().numpy()

y_true = (
    data.y[data.val_mask]
    .cpu()
    .numpy()
)

y_pred = (
    pred[data.val_mask]
    .cpu()
    .numpy()
)

y_prob = (
    prob[data.val_mask]
    .cpu()
    .numpy()
)

print("Confusion matrix:")
print(
    confusion_matrix(
        y_true,
        y_pred
    )
)

# Confusion matrix:
# [[4413    3]
#  [  14  329]]


print("")
print("Classification report:")
print(
    classification_report(
        y_true,
        y_pred
    )
)

# 模型的表現:
# Confusion matrix:
# [[4412    4]
#  [  10  333]]

# Classification report:
#               precision    recall  f1-score   support

#            0       1.00      1.00      1.00      4416
#            1       0.99      0.97      0.98       343

#     accuracy                           1.00      4759
#    macro avg       0.99      0.98      0.99      4759
# weighted avg       1.00      1.00      1.00      4759

print("----------------------------------------------------------------")
print("ROC的結果:")
from sklearn.metrics import roc_curve, roc_auc_score
import matplotlib.pyplot as plt


# 計算 ROC curve
fpr, tpr, thresholds = roc_curve(
    y_true,
    y_prob
)


# 計算 AUC
auc = roc_auc_score(
    y_true,
    y_prob
)


plt.figure(figsize=(6,6))


plt.plot(
    fpr,
    tpr,
    label=f"GCN (AUC={auc:.3f})"
)


# random classifier
plt.plot(
    [0,1],
    [0,1],
    linestyle="--",
    label="Random"
)


plt.xlabel(
    "False Positive Rate"
)

plt.ylabel(
    "True Positive Rate"
)


plt.title(
    "ROC Curve"
)


plt.legend()

plt.show()

print("")

print(
    "ROC-AUC:",
    round(auc,4)
)

# ROC-AUC: 0.9999

from sklearn.metrics import average_precision_score


pr_auc = average_precision_score(
    y_true,
    y_prob
)


print(
    "PR-AUC:",
    round(pr_auc,4)
)

# PR-AUC: 0.9983

# 探討模型給每個細胞的預測結果以及機率
healthy = metadata[metadata["label"] == 0]

infected = metadata[metadata["label"] == 1]

unknown_pred_healthy = metadata[
    (metadata["label"] == -1) &
    (metadata["predict_label"] == 0)
]

unknown_pred_infected = metadata[
    (metadata["label"] == -1) &
    (metadata["predict_label"] == 1)
]

unknown = metadata[metadata["label"] == -1]

# 將原本未知的細胞資訊匯出
# save_path = r"C:\Users\USER\Desktop\資訊所實習\計畫\資料探勘\Ravindra2021.raw_count.stdprep.h5ad的細胞標籤\GCN_unknown_prediction.csv"

# unknown.to_csv(
#     save_path,
#     encoding="utf-8-sig"
# )

# print("CSV 已輸出:")
# print(save_path)
# print("資料大小:", unknown.shape)

plt.figure(figsize=(8,6))

plt.hist(
    np.log10(infected["Viral_transcript"] + 1),
    bins=50,
    alpha=0.5,
    label="Original infected"
)

plt.hist(
    np.log10(unknown_pred_healthy["Viral_transcript"] + 1),
    bins=50,
    alpha=0.5,
    label="Unknown -> Healthy"
)

plt.hist(
    np.log10(unknown_pred_infected["Viral_transcript"] + 1),
    bins=50,
    alpha=0.5,
    label="Unknown -> Infected"
)

plt.axvline(
    np.log10(11),   # 對應 Viral transcript = 10
    color="red",
    linestyle="--",
    linewidth=2,
    label="Paper cutoff = 10"
)

plt.yscale("log")   # y 軸改成 log scale

plt.xlabel("log10(Viral transcript + 1)")
plt.ylabel("Cell count (log scale)")

plt.legend()
plt.tight_layout()
plt.show()

print("")
print("----------------------------------------------------------------")
print("未知的細胞最後被判定為感染的summary:")

low_viral = unknown_pred_infected[
    unknown_pred_infected["Viral_transcript"] < 10
]

print(f"Unknown → Infected: {len(unknown_pred_infected)} cells")
print(f"Viral transcript < 10: {len(low_viral)} cells")
print(f"Percentage: {len(low_viral)/len(unknown_pred_infected)*100:.2f}%")
# 未知的細胞最後被判定為感染的summary:
# Unknown → Infected: 37282 cells
# Viral transcript < 10: 34911 cells
# Percentage: 93.64%


print("病毒轉錄體的分布:")
viral = unknown_pred_infected["Viral_transcript"]

print("0:", (viral == 0).sum())
print("1~9:", ((viral > 0) & (viral < 10)).sum())
print("10~99:", ((viral >= 10) & (viral < 100)).sum())
print(">=100:", (viral >= 100).sum())

# 病毒轉錄體的分布:
# 0: 8947
# 1~9: 25964
# 10~99: 1906
# >=100: 465


print("")
print("模型判斷的信心:")
print(unknown_pred_infected["infect_prob"].describe())

# 模型判斷的信心:
# count    37282.000000
# mean         0.929394
# std          0.123776
# min          0.500053
# 25%          0.918473
# 50%          0.997383
# 75%          0.999997
# max          1.000000

print("")
print("根本沒有病毒轉錄體，但被判定為感染的信心:")
zero_cells = unknown_pred_infected[
    unknown_pred_infected["Viral_transcript"] == 0
]

print(
    zero_cells["infect_prob"].describe()
)

# 根本沒有病毒轉錄體，但被判定為感染的信心:
# count    8947.000000
# mean        0.920629
# std         0.129182
# min         0.500053
# 25%         0.892943
# 50%         0.994462
# 75%         0.999963
# max         1.000000

# 檢查病毒轉錄體是否仍然是判斷感染與否的標準
import seaborn as sns
import matplotlib.pyplot as plt

plot_df = metadata.copy()

plot_df["Transcript_group"] = pd.cut(
    plot_df["Viral_transcript"],
    bins=[-1,0,9,99,float("inf")],
    labels=["0","1-9","10-99",">=100"]
)

plt.figure(figsize=(8,5))

sns.boxplot(
    data=plot_df,
    x="Transcript_group",
    y="infect_prob",
    showfliers=False
)

plt.xlabel("Viral transcript")
plt.ylabel("Predicted infection probability")

plt.show()