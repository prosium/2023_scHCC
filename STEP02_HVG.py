import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import scanpy as sc
import decoupler as dc
import numpy as np
import scipy.sparse as sp

all_data.obs_names_make_unique()
all_data

if sp.issparse(all_data.X):
    all_data.X = all_data.X.toarray()

all_data.X = np.nan_to_num(all_data.X)
all_data.X[all_data.X == np.inf] = 0
all_data.X[all_data.X == -np.inf] = 0
all_data.X[all_data.X > 1e6] = 0

sc.pp.normalize_total(all_data, target_sum=1e4)
sc.pp.log1p(all_data)


all_data.raw = all_data
sc.pp.highly_variable_genes(all_data, min_mean=0.0125, max_mean=3, min_disp=0.5)
all_data = all_data[:, all_data.var['highly_variable']]

print("Highly variable genes: %d"%sum(all_data.var.highly_variable))

#plot variable genes
sc.pl.highly_variable_genes(all_data)

# subset for variable genes in the dataset
all_data = all_data[:, all_data.var['highly_variable']]
print(all_data.X.shape)
print(all_data.raw.X.shape)

sc.pp.regress_out(all_data, ['total_counts', 'pct_counts_mt','pct_counts_ribo'])
sc.pp.scale(all_data, max_value=10)

sc.tl.pca(all_data, svd_solver = 'arpack')
