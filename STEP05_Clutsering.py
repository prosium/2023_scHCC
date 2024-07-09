import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import decoupler as dc

sc.tl.leiden(all_data, key_added = "leiden_1.0", resolution=1.0)

sc.settings.set_figure_params(dpi=200, figsize=(4,4))
sc.pl.tsne(all_data, color='leiden_1.0')

dc.run_ora(
    mat=all_data,
    net=markers,
    source='cell_type',
    target='genesymbol',
    min_n=3,
    verbose=True,
    use_raw=False
)

acts = dc.get_acts(all_data, obsm_key='ora_estimate')

# We need to remove inf and set them to the maximum value observed for pvals=0
acts_v = acts.X.ravel()
max_e = np.nanmax(acts_v[np.isfinite(acts_v)])
acts.X[~np.isfinite(acts.X)] = max_e

df = dc.rank_sources_groups(acts, groupby='leiden_1.0', reference='rest', method='t-test_overestim_var')
n_ctypes = 3
ctypes_dict = df.groupby('group').head(n_ctypes).groupby('group')['names'].apply(lambda x: list(x)).to_dict()
annotation_dict = df.groupby('group').head(1).set_index('group')['names'].to_dict()
all_data.obs['Azimuth_L2'] = [annotation_dict[clust] for clust in all_data.obs['leiden_1.0']]

sc.pl.tsne(all_data, color='Azimuth_L2')

import seaborn as sns
import matplotlib.pyplot as plt
cell_counts = all_data.obs.groupby('Azimuth_L1').size().reset_index(name='cell_count')


print(cell_counts)

cluster_counts = all_data.obs.groupby('Azimuth_L1')['leiden_1.0'].nunique().reset_index(name='n_clusters')


print(cluster_counts)

# subset_endothelials = all_data[all_data.obs['Azimuth_L1'] == 'Endothelials']
marker = pd.read_table("../Marker_Azimuth_L2.txt", sep="\t")
cell_types_of_interest = ["Sinusoidal Endothelial", "Vascual Central Venous Endothelial", "Vascual Portal Endothelial"]

filtered_marker = marker[marker['cell_type'].isin(cell_types_of_interest)]








