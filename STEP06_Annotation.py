import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import decoupler as dc

all_data.obs['Azimuth_L1'].value_counts()
cell_types_of_interest = ['Hepatic stellate cells']
filtered_data = all_data[all_data.obs['Azimuth_L1'].isin(cell_types_of_interest), :]
print(filtered_data)

# subset_endothelials = all_data[all_data.obs['Azimuth_L1'] == 'Endothelials']
marker = pd.read_table("../Marker_Azimuth_L2.txt", sep="\t")
cell_types_of_interest = ["Hepatic stellate cells", "Myofibroblast", "Portal Fibroblast"]
filtered_marker = marker[marker['cell_type'].isin(cell_types_of_interest)]
print(filtered_marker['cell_type'].value_counts())

dc.run_ora(
    mat=filtered_data,
    net=filtered_marker,
    source='cell_type',
    target='genesymbol',
    verbose=True,
    use_raw=False
)

acts = dc.get_acts(filtered_data, obsm_key='ora_estimate')

# We need to remove inf and set them to the maximum value observed for pvals=0
acts_v = acts.X.ravel()
max_e = np.nanmax(acts_v[np.isfinite(acts_v)])
acts.X[~np.isfinite(acts.X)] = max_e

df = dc.rank_sources_groups(acts, groupby='leiden_1.0', reference='rest', method='t-test_overestim_var')
n_ctypes = 3
ctypes_dict = df.groupby('group').head(n_ctypes).groupby('group')['names'].apply(lambda x: list(x)).to_dict()
annotation_dict = df.groupby('group').head(1).set_index('group')['names'].to_dict()
filtered_data.obs['Azimuth_L2'] = [annotation_dict[clust] for clust in filtered_data.obs['leiden_1.0']]

sc.pl.tsne(filtered_data, color="Azimuth_L2")

all_data.obs.loc[filtered_data.obs.index, 'Azimuth_L2'] = filtered_data.obs['Azimuth_L2']
# all_data.obs.loc[all_data.obs.index, 'Azimuth_L2'] = all_data.obs['Azimuth_L1']
all_data.obs.loc[all_data.obs['Azimuth_L1'] == 'B cells', 'Azimuth_L2'] = 'B cells'


