import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import warnings
import os
import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")


obs_df = all_data.obs
obs_df['Group'] = pd.Categorical(obs_df['Group'], categories=['Normal', 'Hepatitis', 'HCC'], ordered=True)

grouped = obs_df.groupby(['Group', 'Azimuth_L2']).size().unstack(fill_value=0)


proportions = grouped.div(grouped.sum(axis=1), axis=0)

extended_palette = ['#1F77B4', '#FF7F0E', '#2CA02C', '#D62728', '#9467BD',
                    '#8C564B', '#E377C2', '#7F7F7F', '#BCBD22', '#17BECF',
                    '#AEC7E8', '#FFBB78', '#98DF8A', '#FF9896']

ax = proportions.plot(kind='bar', stacked=True, figsize=(7, 7), color=extended_palette, width=0.8)

plt.title('Stacked Bar Plot of Azimuth_L2 Proportions by From')
plt.xlabel('From')
plt.ylabel('Proportion')
plt.xticks(rotation=90)
plt.legend(title='Azimuth_L2', bbox_to_anchor=(1.05, 1), loc='upper left')
ax.xaxis.grid(False)
plt.tight_layout()
plt.show()


import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd

obs_df = all_data.obs
grouped = obs_df.groupby(['Azimuth_L2', 'from']).size().unstack(fill_value=0)
extended_palette = ['#1F77B4', '#FF7F0E', '#2CA02C', '#D62728', '#9467BD',
                    '#8C564B', '#E377C2', '#7F7F7F', '#BCBD22', '#17BECF',
                    '#AEC7E8', '#FFBB78', '#98DF8A', '#FF9896']

fig, ax = plt.subplots(figsize=(7, 7), dpi=300)
grouped.plot(kind='bar', stacked=True, color=extended_palette, ax=ax, width=0.8)
plt.title('Stacked Bar Plot of Azimuth_L2 Counts by From')
plt.xlabel('From')
plt.ylabel('Count')
plt.xticks(rotation=90)
plt.legend(title='Azimuth_L2', bbox_to_anchor=(1.05, 1), loc='upper left')

ax.xaxis.grid(False)

plt.tight_layout()
plt.show()


import pandas as pd
import matplotlib.pyplot as plt

obs_df = all_data.obs
obs_df['Group'] = pd.Categorical(obs_df['Tissue'], categories=['Normal', 'Tumor'], ordered=True)
grouped = obs_df.groupby(['Tissue', 'Azimuth_L2']).size().unstack(fill_value=0)
proportions = grouped.div(grouped.sum(axis=1), axis=0)
extended_palette = ['#1F77B4', '#FF7F0E', '#2CA02C', '#D62728', '#9467BD',
                    '#8C564B', '#E377C2', '#7F7F7F', '#BCBD22', '#17BECF',
                    '#AEC7E8', '#FFBB78', '#98DF8A', '#FF9896']

fig, ax = plt.subplots(figsize=(7.5, 7))


for idx, azimuth in enumerate(proportions.columns):
    ax.plot(proportions.index, proportions[azimuth], marker='o', color=extended_palette[idx % len(extended_palette)], label=azimuth)
plt.title('Line Plot of Azimuth_L2 Proportions by Group')
plt.xlabel('Group')
plt.ylabel('Proportion')
plt.xticks(rotation=45)
plt.legend(title='Azimuth_L2', bbox_to_anchor=(1.05, 1), loc='upper left')
ax.xaxis.grid(False)

plt.tight_layout()
plt.show()
