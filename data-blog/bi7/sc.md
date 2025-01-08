# Density based clustering parametric test

## Load package

```python
import pandas as pd
import numpy as np
import os
import math
import matplotlib.pyplot as plt
```

## Set path

```python
os.chdir("/data-blog/bi7")
os.getcwd()
```
```plain text
'/data-blog/bi7'
```

## Parametric test

```python
k = 5

# Load input data
input_df = pd.read_csv('mutClustInput_df.csv')
cluster477_df = pd.read_csv('clusters_min-clust-5.tsv', sep='\t')
cluster28_df = pd.read_csv('Clusters_28.csv')

# Add mutation_count column
cluster477_df['mutation_count'] = cluster477_df['mut_positions'].str.split(',').apply(len)

# Add eps column
cluster477_df['eps'] = cluster477_df['length'].apply(lambda x: math.floor((x - 1) / 2))

# Filter input data within a specific range
filtered_input_df = input_df[(input_df.Position >= 266) & (input_df.Position <= 29674)]

# Calculate cumulative sum and count for 'H-score'
input_df['cumsum_hscore'] = input_df['H-score'].cumsum()
input_df['cumcount_hscore'] = input_df['H-score'].ne(0).cumsum()

# Define a function to get sum and count within a range
def get_sum_count(pos, eps):
    lower_bound = max(0, pos - eps - 266)
    upper_bound = min(pos + eps - 266, len(input_df) - 1)
    hsum = input_df['cumsum_hscore'].iloc[upper_bound] - input_df['cumsum_hscore'].iloc[lower_bound]
    hcount = input_df['cumcount_hscore'].iloc[upper_bound] - input_df['cumcount_hscore'].iloc[lower_bound]
    return hsum, hcount

# Calculate eps values and lengths
list_eps = []
list_length = []
for pos in filtered_input_df['Position']:
    eps = 1
    hmean, hsum, mtcount = 0, 0, 0

    while hmean < 0.01 or hsum < 0.05 or mtcount < k:
        if eps > 1000:
            eps = 1000
            break

        hsum, mtcount = get_sum_count(pos, eps)
        hmean = hsum / ((2 * eps) + 1)
        eps += 1

    list_eps.append(eps)
    list_length.append(min(eps * 2 + 1, 2001))

filtered_input_df[f'{k}-dist'] = list_length

# Determine cluster membership for cluster_477
filtered_input_df['cluster_477'] = filtered_input_df['Position'].apply(
    lambda pos: 'Y' if any((pos >= left) and (pos <= right) for left, right in 
                           zip(cluster477_df['left_position'], cluster477_df['right_position'])) else 'N'
)

# Determine cluster membership for cluster_28
filtered_input_df['cluster_28'] = filtered_input_df['Position'].apply(
    lambda pos: 'Y' if any((pos >= left) and (pos <= right) for left, right in 
                           zip(cluster28_df['left_position'], cluster28_df['right_position'])) else 'N'
)

# Save and sort filtered data
filtered_input_df = filtered_input_df.sort_values(by=f'{k}-dist', ascending=False).reset_index(drop=True)
#filtered_input_df.to_csv(f"input_df_{k}-dist.tsv", sep='\t', index=False)
```
```python
filtered_input_df['cluster_477'].value_counts()
filtered_input_df['cluster_28'].value_counts()
```
```
N    21007
Y     8402
Name: cluster_477, dtype: int64

N    28168
Y     1241
Name: cluster_28, dtype: int64
```

### Visualization

```python
for index, row in filtered_input_df.iterrows():
    if row['cluster_477'] == 'Y':
        if row['cluster_28'] == 'Y':
            plt.scatter(index, row[f'{k}-dist'], s=3, c='yellow') # 28인 경우 yellow
        else:
            plt.scatter(index, row[f'{k}-dist'], s=1, c='blue')  # 477인 경우 blue
    else:
        plt.scatter(index, row[f'{k}-dist'], s=1, c='lightgray')  # N인 경우 lightgray
plt.xlabel('Index')
plt.ylabel(f'{k}-dist')
plt.title(f'{k}-dist plot')
plt.grid(True)
plt.show()
```
![image](https://github.com/user-attachments/assets/2d653fcf-f0cb-4f8f-b382-687ba7c7fb09)
