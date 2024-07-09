import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

all_data.var_names
all_data.uns['log1p']['base'] = None
sc.pp.highly_variable_genes(all_data, min_mean = 0.0125, max_mean = 3, min_disp = 0.5, batch_key = 'from')

var_genes_batch = all_data.var.highly_variable_nbatches > 0
var_select = all_data.var.highly_variable_nbatches > 1
var_genes = var_select.index[var_select]
batches = all_data.obs['from'].cat.categories.tolist()

Alldata = {}
for batch in batches:
        Alldata[batch] = all_data[all_data.obs['from'] == batch,]

import scanorama
Alldata2 = dict()
for ds in Alldata.keys():
        print(ds)
        Alldata2[ds] = Alldata[ds][:,var_genes]
adatas = list(Alldata2.values())
scanorama.integrate_scanpy(adatas, dimred = 50)

scanorama_int = [ad.obsm['X_scanorama'] for ad in adatas]
all_s = np.concatenate(scanorama_int)
adata_sc = all_data.copy()
adata_sc.obsm["Scanorama"] = all_s

adata_sc.write('batchcorrected.h5ad')
