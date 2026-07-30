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

# Classification report:
#               precision    recall  f1-score   support

#            0       1.00      1.00      1.00      4416
#            1       0.99      0.96      0.97       343

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

# PR-AUC: 0.9986