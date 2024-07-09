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

all_data.obs['Azimuth_L2'].value_counts
sc.tl.rank_genes_groups(all_data, groupby="Azimuth_L2", method="t-test_overestim_var")
sc.pl.rank_genes_groups_dotplot(
    all_data, groupby="Azimuth_L2", standard_scale="var", n_genes=5)
