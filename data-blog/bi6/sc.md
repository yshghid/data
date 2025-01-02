# Random forest GC subtype 분류기 생성

## Load package

```python
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns

import sklearn
from sklearn.preprocessing import QuantileTransformer
from sklearn.datasets import make_blobs
from sklearn.manifold import TSNE
from scipy.stats import nbinom

import joypy
import random
from matplotlib import cm

from sklearn.preprocessing import quantile_transform
from sklearn.preprocessing import MinMaxScaler

from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, classification_report
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix
```

## Set path

```python
os.chdir("/data-blog/bi6")
os.getcwd()
```
```plain text
'/data-blog/bi6'
```

## RF Classifier

Normalize data

```python
count = pd.read_csv("count.csv", encoding='cp949')
count_mat = count.set_index('gene_symbol').drop(columns='ensembl_id')
count_log = np.log2(count_mat + 1)
count_T = count_mat.T

meta = pd.read_csv("meta.csv",encoding='cp949')
meta['Subtype'] = meta['Subtype'].apply(lambda x: 'other' if x != 'EBV' else x)
meta_mat = meta.set_index('sampleID')
meta_mat = meta_mat.drop('Hypermethylation-category', axis=1)
meta_mat['Subtype'].value_counts()
subtypes = meta_mat['Subtype'].tolist()
```

Draw TSNE: raw data vs log2 normalized data

```python
# tsne 1
tsne = TSNE(n_components=2, perplexity=30, random_state=42)
data_tsne = tsne.fit_transform(count_T)
tsne_df = pd.DataFrame(data_tsne, columns=['x', 'y'])
tsne_df['Subtype'] = subtypes
plt.figure(figsize=(10, 8))
subtype_colors = {'other': 'red', 'EBV': 'green'}
subtype_labels = list(subtype_colors.keys())
for subtype in subtype_labels:
    subtype_indices = tsne_df[tsne_df['Subtype'] == subtype]
    plt.scatter(subtype_indices['x'], subtype_indices['y'], color=subtype_colors[subtype], label=subtype)
plt.title('t-SNE Plot')
plt.xlabel('t-SNE Component 1')
plt.ylabel('t-SNE Component 2')
plt.legend(title='Subtypes')
plt.show()

# tsne 2
subtypes = meta_mat['Subtype'].tolist()
count_log_T = count_log.T
tsne = TSNE(n_components=2, perplexity=30, random_state=42)
data_tsne = tsne.fit_transform(count_log_T)
tsne_df = pd.DataFrame(data_tsne, columns=['x', 'y'])
tsne_df['Subtype'] = subtypes
plt.figure(figsize=(10, 8))
subtype_colors = {'other': 'red', 'EBV': 'green'}
subtype_labels = list(subtype_colors.keys())
for subtype in subtype_labels:
    subtype_indices = tsne_df[tsne_df['Subtype'] == subtype]
    plt.scatter(subtype_indices['x'], subtype_indices['y'], color=subtype_colors[subtype], label=subtype)
plt.title('t-SNE Plot')
plt.xlabel('t-SNE Component 1')
plt.ylabel('t-SNE Component 2')
plt.legend(title='Subtypes')
plt.show()
```

Draw distribution: raw data vs log2 transformed data

```python
count_log_random = count_log.iloc[random.sample(range(len(count_log)), 1000)]
count_random = count_mat.iloc[random.sample(range(len(count_mat)), 1000)]

# hist1
sns.histplot(count_random.iloc[:,0], kde=True)
plt.title("Distribution of data")
plt.xlabel("Value")
plt.ylabel("Frequency")
plt.show()

# hist2
sns.histplot(count_log_random.iloc[:,0], kde=True)
plt.title("Distribution of log2 transformed data")
plt.xlabel("Value")
plt.ylabel("Frequency")
plt.show()
```

Random sampling

```python
count_log_T = count_log.T
count_log_T['sampleID'] = count_log_T.index

meta2 = meta.drop("Hypermethylation-category", axis=1)

data = pd.merge(count_log_T, meta2, on='sampleID')

# random sampling 40 samples
other_rows = data[data['Subtype'] == 'other']
random_indices = random.sample(list(other_rows.index), 40)
data = pd.concat([data.loc[random_indices], data.loc[data['Subtype'] == 'EBV']])

meta2['Subtype'].value_counts()
data['Subtype'].value_counts()
```
```plain text
other    315
EBV       30
Name: Subtype, dtype: int64

other    40
EBV      30
Name: Subtype, dtype: int64
```

RF classifier

```python
X = data.iloc[:, :-1]  # gene expression 데이터
y = data['Subtype']   # subtype을 예측할 label

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
rf_classifier = RandomForestClassifier(n_estimators=100, random_state=42)#, max_depth=20)
rf_classifier.fit(X_train.iloc[:, :-2], y_train)
y_pred = rf_classifier.predict(X_test.iloc[:, :-2])
accuracy = accuracy_score(y_test, y_pred)

print("Accuracy:", accuracy)
print("Classification Report:")
print(classification_report(y_test, y_pred))

conf_matrix = confusion_matrix(y_test, y_pred)
confusion_matrix_df = pd.DataFrame(conf_matrix, index=rf_classifier.classes_, columns=rf_classifier.classes_)
print(confusion_matrix_df)
```
```plain text
Accuracy: 0.8571428571428571
Classification Report:
              precision    recall  f1-score   support

         EBV       0.67      1.00      0.80         4
       other       1.00      0.80      0.89        10

    accuracy                           0.86        14
   macro avg       0.83      0.90      0.84        14
weighted avg       0.90      0.86      0.86        14

       EBV  other
EBV      4      0
other    2      8
```

Extract feature importance

```python
feature_importance = rf_classifier.feature_importances_
importance_df = pd.DataFrame({'Feature': X_train.columns[:-2], 'Importance': feature_importance})
importance_df = importance_df.sort_values(by='Importance', ascending=False)
importance_df_head = importance_df.head(100)
list_head = importance_df_head['Feature'].tolist()
print(*list_head)
```
```plain text
ERVMER34-1 FCGR1A COX6B2 ITGAE WDR72 IKZF3 PRAC2 FCGR3A CAMK2N2 ADAMDEC1 APOBEC3H CAMK2N1 CLEC12A CLDN3 NOL3 NKG7 HYLS1 CRAT KCNH8 CXCL11 GZMA DUOX2 TXNRD3 GBP4 LGMN SSR4 GNGT2 CD8A KLHL35 FNDC11 PPIE HLA-DMB PDIA2 GZMB TOR4A FAM110A IPCEF1 SYT1 TMEM86A HOXA10 DDB2 TRIAP1 SLA2 CCNP YBX2 PRKCQ NYAP1 BECN2 ZNF683 FAH RAD51C GOLT1A KLRK1 RGS1 ENTPD8 PLAG1 DUOXA2 ASS1 SLC5A12 BCL2L10 HEBP1 CTSC TDRD5 GJB5 OR2I1P PAQR7 HLA-DPA1 FTL ARHGEF3 ADPRHL1 NAB1 PTGS2 FAIM CXCR6 CHAF1B CABYR DRD1 HOXC4 LRRIQ1 ITGB7 HOXA5 PDE4A RNF208 C19orf38 MBD4 TNNC1 WNT16 PRXL2A ADGRG1 KCNK15 ZNF491 SLC6A20 HOXB6 SLC44A5 NT5DC3 TRIM43 FXYD5 PLPP1 ABI3 TP53I13
```

Draw TSNE: Feature extracted data
```python
count_log_T_head = count_log_T[list_head]
count_log_T_head_T = count_log_T_head.T

# tsne 3
subtypes = meta_mat['Subtype'].tolist()

tsne = TSNE(n_components=2, perplexity=30, random_state=42)
data_tsne = tsne.fit_transform(count_log_T_head)

tsne_df = pd.DataFrame(data_tsne, columns=['x', 'y'])
tsne_df['Subtype'] = subtypes

plt.figure(figsize=(10, 8))
subtype_colors = {'other': 'red', 'EBV': 'green'}
subtype_labels = list(subtype_colors.keys())
for subtype in subtype_labels:
    subtype_indices = tsne_df[tsne_df['Subtype'] == subtype]
    plt.scatter(subtype_indices['x'], subtype_indices['y'], color=subtype_colors[subtype], label=subtype)
plt.title('t-SNE Plot')
plt.xlabel('t-SNE Component 1')
plt.ylabel('t-SNE Component 2')
plt.legend(title='Subtypes')
plt.show()
```
