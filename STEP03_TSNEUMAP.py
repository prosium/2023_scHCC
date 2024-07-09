import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

sc.tl.tsne(all_data)
sc.pp.neighbors(all_data, n_pcs = 30, n_neighbors = 20)
sc.tl.umap(all_data)

order = ['Ours', 'CRA002308', 'GSE156625', 'GSE182159', 'GSE186343', 'GSE234241', 'GSE252863', 'SRP318499']
all_data.obs['from'] = pd.Categorical(all_data.obs['from'], categories=order, ordered=True)
sc.pl.tsne(all_data, color='from')

sc.pl.umap(all_data, color='from')

import scanpy as sc
import matplotlib.pyplot as plt

from_values = all_data.obs['from'].unique()
n_cols = 2
n_rows = (len(from_values) + n_cols - 1) // n_cols  # 총 행 수 계산

fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 5 * n_rows))

for i, from_value in enumerate(from_values):
    ax = axes[i // n_cols, i % n_cols]
    subset = all_data[all_data.obs['from'] == from_value]
    sc.pl.umap(subset, title=f'UMAP for {from_value}', show=False, ax=ax)


plt.tight_layout()
plt.show()

cell_counts = all_data.obs['from'].value_counts()

cell_counts_df = pd.DataFrame(cell_counts).reset_index()
cell_counts_df.columns = ['from', 'cell_count']

all_data2 = all_data.copy()
all_data2.uns['log1p']['base'] = None
sc.pp.highly_variable_genes(all_data2, min_mean = 0.0125, max_mean = 3, min_disp = 0.5, batch_key = 'from')

print("Highly variable genes intersection: %d"%sum(all_data2.var.highly_variable_intersection))
print(all_data2.var.highly_variable_nbatches.value_counts())
var_genes_batch = all_data2.var.highly_variable_nbatches > 0
var_genes_batch
var_select = all_data2.var.highly_variable_nbatches > 1
var_genes = var_select.index[var_select]
batches = all_data.obs['from'].cat.categories.tolist()

Alldata = {}
for batch in batches:
        Alldata[batch] = all_data2[all_data2.obs['from'] == batch,]

Alldata
Alldata2 = dict()
for ds in Alldata.keys():
        print(ds)
        Alldata2[ds] = Alldata[ds][:,var_genes]
adatas = list(Alldata2.values())
scanorama.integrate_scanpy(adatas, dimred = 50)

