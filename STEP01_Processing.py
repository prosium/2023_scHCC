import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import warnings
import os
warnings.filterwarnings("ignore")

SOLD001=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD001.h5ad")
SOLD002=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD002.h5ad")
SOLD003=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD003.h5ad")
SOLD004=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD004.h5ad")
SOLD005=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD005.h5ad")
SOLD006=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD006.h5ad")
SOLD007=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD007.h5ad")
SOLD008=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD008.h5ad")
SOLD009=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD009.h5ad")
SOLD010=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD010.h5ad")
SOLD011=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD011.h5ad")
SOLD012=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD012.h5ad")
SOLD013=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD013.h5ad")
SOLD014=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD014.h5ad")
SOLD015=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD015.h5ad")
SOLD016=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD016.h5ad")
SOLD017=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD017.h5ad")
SOLD018=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD018.h5ad")
SOLD019=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD019.h5ad")
SOLD020=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD020.h5ad")
SOLD021=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD021.h5ad")
SOLD022=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD022.h5ad")
SOLD023=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD023.h5ad")
SOLD024=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD024.h5ad")
SOLD025=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD025.h5ad")
SOLD026=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD026.h5ad")
SOLD027=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD027.h5ad")
SOLD028=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD028.h5ad")
SOLD029=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD029.h5ad")
SOLD030=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD030.h5ad")
SOLD031=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD031.h5ad")
SOLD032=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD032.h5ad")
SOLD033=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD033.h5ad")
SOLD034=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD034.h5ad")
SOLD035=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD035.h5ad")
SOLD036=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD036.h5ad")
SOLD037=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD037.h5ad")
SOLD038=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD038.h5ad")
SOLD039=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD039.h5ad")
SOLD040=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD040.h5ad")
SOLD041=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD041.h5ad")
SOLD042=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD042.h5ad")
SOLD043=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD043.h5ad")
SOLD044=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD044.h5ad")
SOLD045=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD045.h5ad")
SOLD046=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD046.h5ad")
SOLD047=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD047.h5ad")
SOLD048=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD048.h5ad")
SOLD049=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD049.h5ad")
SOLD050=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD050.h5ad")
SOLD051=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD051.h5ad")
SOLD052=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD052.h5ad")
SOLD053=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD053.h5ad")
SOLD054=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD054.h5ad")
SOLD055=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD055.h5ad")
SOLD056=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD056.h5ad")
SOLD057=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD057.h5ad")
SOLD058=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD058.h5ad")
SOLD059=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD059.h5ad")
SOLD060=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD060.h5ad")
SOLD061=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD061.h5ad")
SOLD062=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD062.h5ad")
SOLD063=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD063.h5ad")
SOLD064=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD064.h5ad")
SOLD065=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD065.h5ad")
SOLD066=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD066.h5ad")
SOLD067=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD067.h5ad")
SOLD068=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD068.h5ad")
SOLD069=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD069.h5ad")
SOLD070=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD070.h5ad")
SOLD071=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD071.h5ad")
SOLD072=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD072.h5ad")
SOLD073=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD073.h5ad")
SOLD074=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD074.h5ad")
SOLD075=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD075.h5ad")
SOLD076=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD076.h5ad")
SOLD077=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD077.h5ad")
SOLD078=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD078.h5ad")
SOLD079=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD079.h5ad")
SOLD080=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD080.h5ad")
SOLD081=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD081.h5ad")
SOLD082=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD082.h5ad")
SOLD083=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD083.h5ad")
SOLD084=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD084.h5ad")
SOLD085=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD085.h5ad")
SOLD086=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD086.h5ad")
SOLD087=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD087.h5ad")
SOLD088=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD088.h5ad")
SOLD089=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD089.h5ad")
SOLD090=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD090.h5ad")
SOLD091=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD091.h5ad")
SOLD092=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD092.h5ad")
SOLD093=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD093.h5ad")
SOLD094=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD094.h5ad")
SOLD095=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD095.h5ad")
SOLD096=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD096.h5ad")
SOLD097=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD097.h5ad")
SOLD098=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD098.h5ad")
SOLD099=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD099.h5ad")
SOLD100=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD100.h5ad")
SOLD101=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD101.h5ad")
SOLD102=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD102.h5ad")
SOLD103=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD103.h5ad")
SOLD104=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD104.h5ad")
SOLD105=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD105.h5ad")
SOLD106=sc.read_h5ad("../../SCDATA/01.Datasets/SOLD106.h5ad")


SOLD001.var_names_make_unique()
SOLD001.var_names_make_unique()
SOLD002.var_names_make_unique()
SOLD002.var_names_make_unique()
SOLD003.var_names_make_unique()
SOLD003.var_names_make_unique()
SOLD004.var_names_make_unique()
SOLD004.var_names_make_unique()
SOLD005.var_names_make_unique()
SOLD005.var_names_make_unique()
SOLD006.var_names_make_unique()
SOLD006.var_names_make_unique()
SOLD007.var_names_make_unique()
SOLD007.var_names_make_unique()
SOLD008.var_names_make_unique()
SOLD008.var_names_make_unique()
SOLD009.var_names_make_unique()
SOLD009.var_names_make_unique()
SOLD010.var_names_make_unique()
SOLD010.var_names_make_unique()
SOLD011.var_names_make_unique()
SOLD011.var_names_make_unique()
SOLD012.var_names_make_unique()
SOLD012.var_names_make_unique()
SOLD013.var_names_make_unique()
SOLD013.var_names_make_unique()
SOLD014.var_names_make_unique()
SOLD014.var_names_make_unique()
SOLD015.var_names_make_unique()
SOLD015.var_names_make_unique()
SOLD016.var_names_make_unique()
SOLD016.var_names_make_unique()
SOLD017.var_names_make_unique()
SOLD017.var_names_make_unique()
SOLD018.var_names_make_unique()
SOLD018.var_names_make_unique()
SOLD019.var_names_make_unique()
SOLD019.var_names_make_unique()
SOLD020.var_names_make_unique()
SOLD020.var_names_make_unique()
SOLD021.var_names_make_unique()
SOLD021.var_names_make_unique()
SOLD022.var_names_make_unique()
SOLD022.var_names_make_unique()
SOLD023.var_names_make_unique()
SOLD023.var_names_make_unique()
SOLD024.var_names_make_unique()
SOLD024.var_names_make_unique()
SOLD025.var_names_make_unique()
SOLD025.var_names_make_unique()
SOLD026.var_names_make_unique()
SOLD026.var_names_make_unique()
SOLD027.var_names_make_unique()
SOLD027.var_names_make_unique()
SOLD028.var_names_make_unique()
SOLD028.var_names_make_unique()
SOLD029.var_names_make_unique()
SOLD029.var_names_make_unique()
SOLD030.var_names_make_unique()
SOLD030.var_names_make_unique()
SOLD031.var_names_make_unique()
SOLD031.var_names_make_unique()
SOLD032.var_names_make_unique()
SOLD032.var_names_make_unique()
SOLD033.var_names_make_unique()
SOLD033.var_names_make_unique()
SOLD034.var_names_make_unique()
SOLD034.var_names_make_unique()
SOLD035.var_names_make_unique()
SOLD035.var_names_make_unique()
SOLD036.var_names_make_unique()
SOLD036.var_names_make_unique()
SOLD037.var_names_make_unique()
SOLD037.var_names_make_unique()
SOLD038.var_names_make_unique()
SOLD038.var_names_make_unique()
SOLD039.var_names_make_unique()
SOLD039.var_names_make_unique()
SOLD040.var_names_make_unique()
SOLD040.var_names_make_unique()
SOLD041.var_names_make_unique()
SOLD041.var_names_make_unique()
SOLD042.var_names_make_unique()
SOLD042.var_names_make_unique()
SOLD043.var_names_make_unique()
SOLD043.var_names_make_unique()
SOLD044.var_names_make_unique()
SOLD044.var_names_make_unique()
SOLD045.var_names_make_unique()
SOLD045.var_names_make_unique()
SOLD046.var_names_make_unique()
SOLD046.var_names_make_unique()
SOLD047.var_names_make_unique()
SOLD047.var_names_make_unique()
SOLD048.var_names_make_unique()
SOLD048.var_names_make_unique()
SOLD049.var_names_make_unique()
SOLD049.var_names_make_unique()
SOLD050.var_names_make_unique()
SOLD050.var_names_make_unique()
SOLD051.var_names_make_unique()
SOLD051.var_names_make_unique()
SOLD052.var_names_make_unique()
SOLD052.var_names_make_unique()
SOLD053.var_names_make_unique()
SOLD053.var_names_make_unique()
SOLD054.var_names_make_unique()
SOLD054.var_names_make_unique()
SOLD055.var_names_make_unique()
SOLD055.var_names_make_unique()
SOLD056.var_names_make_unique()
SOLD056.var_names_make_unique()
SOLD057.var_names_make_unique()
SOLD057.var_names_make_unique()
SOLD058.var_names_make_unique()
SOLD058.var_names_make_unique()
SOLD059.var_names_make_unique()
SOLD059.var_names_make_unique()
SOLD060.var_names_make_unique()
SOLD060.var_names_make_unique()
SOLD061.var_names_make_unique()
SOLD061.var_names_make_unique()
SOLD062.var_names_make_unique()
SOLD062.var_names_make_unique()
SOLD063.var_names_make_unique()
SOLD063.var_names_make_unique()
SOLD064.var_names_make_unique()
SOLD064.var_names_make_unique()
SOLD065.var_names_make_unique()
SOLD065.var_names_make_unique()
SOLD066.var_names_make_unique()
SOLD066.var_names_make_unique()
SOLD067.var_names_make_unique()
SOLD067.var_names_make_unique()
SOLD068.var_names_make_unique()
SOLD068.var_names_make_unique()
SOLD069.var_names_make_unique()
SOLD069.var_names_make_unique()
SOLD070.var_names_make_unique()
SOLD070.var_names_make_unique()
SOLD071.var_names_make_unique()
SOLD071.var_names_make_unique()
SOLD072.var_names_make_unique()
SOLD072.var_names_make_unique()
SOLD073.var_names_make_unique()
SOLD073.var_names_make_unique()
SOLD074.var_names_make_unique()
SOLD074.var_names_make_unique()
SOLD075.var_names_make_unique()
SOLD075.var_names_make_unique()
SOLD076.var_names_make_unique()
SOLD076.var_names_make_unique()
SOLD077.var_names_make_unique()
SOLD077.var_names_make_unique()
SOLD078.var_names_make_unique()
SOLD078.var_names_make_unique()
SOLD079.var_names_make_unique()
SOLD079.var_names_make_unique()
SOLD080.var_names_make_unique()
SOLD080.var_names_make_unique()
SOLD081.var_names_make_unique()
SOLD081.var_names_make_unique()
SOLD082.var_names_make_unique()
SOLD082.var_names_make_unique()
SOLD083.var_names_make_unique()
SOLD083.var_names_make_unique()
SOLD084.var_names_make_unique()
SOLD084.var_names_make_unique()
SOLD085.var_names_make_unique()
SOLD085.var_names_make_unique()
SOLD086.var_names_make_unique()
SOLD086.var_names_make_unique()
SOLD087.var_names_make_unique()
SOLD087.var_names_make_unique()
SOLD088.var_names_make_unique()
SOLD088.var_names_make_unique()
SOLD089.var_names_make_unique()
SOLD089.var_names_make_unique()
SOLD090.var_names_make_unique()
SOLD090.var_names_make_unique()
SOLD091.var_names_make_unique()
SOLD091.var_names_make_unique()
SOLD092.var_names_make_unique()
SOLD092.var_names_make_unique()
SOLD093.var_names_make_unique()
SOLD093.var_names_make_unique()
SOLD094.var_names_make_unique()
SOLD094.var_names_make_unique()
SOLD095.var_names_make_unique()
SOLD095.var_names_make_unique()
SOLD096.var_names_make_unique()
SOLD096.var_names_make_unique()
SOLD097.var_names_make_unique()
SOLD097.var_names_make_unique()
SOLD098.var_names_make_unique()
SOLD098.var_names_make_unique()
SOLD099.var_names_make_unique()
SOLD099.var_names_make_unique()
SOLD100.var_names_make_unique()
SOLD100.var_names_make_unique()
SOLD101.var_names_make_unique()
SOLD101.var_names_make_unique()
SOLD102.var_names_make_unique()
SOLD102.var_names_make_unique()
SOLD103.var_names_make_unique()
SOLD103.var_names_make_unique()
SOLD104.var_names_make_unique()
SOLD104.var_names_make_unique()
SOLD105.var_names_make_unique()
SOLD105.var_names_make_unique()
SOLD106.var_names_make_unique()

SOLD001.obs['from2']='Ours_HCC'
SOLD002.obs['from2']='Ours_HCC'
SOLD003.obs['from2']='Ours_HCC'
SOLD004.obs['from2']='Ours_HCC'
SOLD005.obs['from2']='Ours_HCC'
SOLD006.obs['from2']='Ours_HCC'
SOLD007.obs['from2']='GSE156625_Normal'
SOLD008.obs['from2']='GSE156625_Normal'
SOLD009.obs['from2']='GSE156625_Normal'
SOLD010.obs['from2']='GSE156625_Normal'
SOLD011.obs['from2']='GSE156625_Normal'
SOLD012.obs['from2']='GSE156625_Normal'
SOLD013.obs['from2']='GSE156625_Normal'
SOLD014.obs['from2']='GSE156625_Normal'
SOLD015.obs['from2']='GSE156625_Normal'
SOLD016.obs['from2']='GSE156625_HCC'
SOLD017.obs['from2']='GSE156625_HCC'
SOLD018.obs['from2']='GSE156625_HCC'
SOLD019.obs['from2']='GSE156625_HCC'
SOLD020.obs['from2']='GSE156625_HCC'
SOLD021.obs['from2']='GSE156625_HCC'
SOLD022.obs['from2']='GSE156625_HCC'
SOLD023.obs['from2']='GSE156625_HCC'
SOLD024.obs['from2']='GSE156625_HCC'
SOLD025.obs['from2']='GSE156625_HCC'
SOLD026.obs['from2']='GSE156625_HCC'
SOLD027.obs['from2']='GSE156625_HCC'
SOLD028.obs['from2']='GSE156625_HCC'
SOLD029.obs['from2']='GSE156625_HCC'
SOLD030.obs['from2']='GSE156625_HCC'
SOLD031.obs['from2']='GSE156625_HCC'
SOLD032.obs['from2']='GSE156625_HCC'
SOLD033.obs['from2']='GSE156625_HCC'
SOLD034.obs['from2']='GSE156625_HCC'
SOLD035.obs['from2']='GSE156625_HCC'
SOLD036.obs['from2']='GSE156625_HCC'
SOLD037.obs['from2']='GSE156625_HCC'
SOLD038.obs['from2']='GSE156625_HCC'
SOLD039.obs['from2']='GSE156625_HCC'
SOLD040.obs['from2']='GSE156625_HCC'
SOLD041.obs['from2']='GSE156625_HCC'
SOLD042.obs['from2']='GSE156625_HCC'
SOLD043.obs['from2']='CRA002308_Normal'
SOLD044.obs['from2']='CRA002308_HCC'
SOLD045.obs['from2']='CRA002308_Normal'
SOLD046.obs['from2']='CRA002308_HCC'
SOLD047.obs['from2']='CRA002308_Normal'
SOLD048.obs['from2']='CRA002308_HCC'
SOLD049.obs['from2']='CRA002308_Normal'
SOLD050.obs['from2']='CRA002308_HCC'
SOLD051.obs['from2']='CRA002308_Normal'
SOLD052.obs['from2']='CRA002308_HCC'
SOLD053.obs['from2']='CRA002308_Normal'
SOLD054.obs['from2']='CRA002308_HCC'
SOLD055.obs['from2']='SRP318499_HCC'
SOLD056.obs['from2']='SRP318499_HCC'
SOLD057.obs['from2']='SRP318499_HCC'
SOLD058.obs['from2']='SRP318499_HCC'
SOLD059.obs['from2']='SRP318499_HCC'
SOLD060.obs['from2']='SRP318499_HCC'
SOLD061.obs['from2']='SRP318499_HCC'
SOLD062.obs['from2']='GSE182159_Normal'
SOLD063.obs['from2']='GSE182159_Normal'
SOLD064.obs['from2']='GSE182159_Normal'
SOLD065.obs['from2']='GSE182159_Normal'
SOLD066.obs['from2']='GSE182159_Normal'
SOLD067.obs['from2']='GSE182159_Normal'
SOLD068.obs['from2']='GSE182159_Hepatitis'
SOLD069.obs['from2']='GSE182159_Hepatitis'
SOLD070.obs['from2']='GSE182159_Hepatitis'
SOLD071.obs['from2']='GSE182159_Hepatitis'
SOLD072.obs['from2']='GSE182159_Hepatitis'
SOLD073.obs['from2']='GSE182159_Hepatitis'
SOLD074.obs['from2']='GSE182159_Hepatitis'
SOLD075.obs['from2']='GSE182159_Hepatitis'
SOLD076.obs['from2']='GSE182159_Hepatitis'
SOLD077.obs['from2']='GSE182159_Hepatitis'
SOLD078.obs['from2']='GSE182159_Hepatitis'
SOLD079.obs['from2']='GSE182159_Hepatitis'
SOLD080.obs['from2']='GSE182159_Hepatitis'
SOLD081.obs['from2']='GSE182159_Hepatitis'
SOLD082.obs['from2']='GSE182159_Hepatitis'
SOLD083.obs['from2']='GSE182159_Hepatitis'
SOLD084.obs['from2']='GSE182159_Hepatitis'
SOLD085.obs['from2']='GSE252863_Hepatitis'
SOLD086.obs['from2']='GSE252863_Hepatitis'
SOLD087.obs['from2']='GSE252863_Hepatitis'
SOLD088.obs['from2']='GSE252863_Hepatitis'
SOLD089.obs['from2']='GSE252863_Hepatitis'
SOLD090.obs['from2']='GSE234241_Hepatitis'
SOLD091.obs['from2']='GSE234241_Hepatitis'
SOLD092.obs['from2']='GSE234241_Hepatitis'
SOLD093.obs['from2']='GSE234241_Hepatitis'
SOLD094.obs['from2']='GSE234241_Hepatitis'
SOLD095.obs['from2']='GSE234241_Hepatitis'
SOLD096.obs['from2']='GSE234241_Hepatitis'
SOLD097.obs['from2']='GSE234241_Hepatitis'
SOLD098.obs['from2']='GSE234241_Normal'
SOLD099.obs['from2']='GSE234241_Normal'
SOLD100.obs['from2']='GSE234241_Normal'
SOLD101.obs['from2']='GSE186343_Hepatitis'
SOLD102.obs['from2']='GSE186343_Hepatitis'
SOLD103.obs['from2']='GSE186343_Hepatitis'
SOLD104.obs['from2']='GSE186343_Hepatitis'
SOLD105.obs['from2']='GSE186343_Hepatitis'
SOLD106.obs['from2']='GSE186343_Hepatitis'

SOLD001.obs['from']='Ours'
SOLD002.obs['from']='Ours'
SOLD003.obs['from']='Ours'
SOLD004.obs['from']='Ours'
SOLD005.obs['from']='Ours'
SOLD006.obs['from']='Ours'
SOLD007.obs['from']='GSE156625'
SOLD008.obs['from']='GSE156625'
SOLD009.obs['from']='GSE156625'
SOLD010.obs['from']='GSE156625'
SOLD011.obs['from']='GSE156625'
SOLD012.obs['from']='GSE156625'
SOLD013.obs['from']='GSE156625'
SOLD014.obs['from']='GSE156625'
SOLD015.obs['from']='GSE156625'
SOLD016.obs['from']='GSE156625'
SOLD017.obs['from']='GSE156625'
SOLD018.obs['from']='GSE156625'
SOLD019.obs['from']='GSE156625'
SOLD020.obs['from']='GSE156625'
SOLD021.obs['from']='GSE156625'
SOLD022.obs['from']='GSE156625'
SOLD023.obs['from']='GSE156625'
SOLD024.obs['from']='GSE156625'
SOLD025.obs['from']='GSE156625'
SOLD026.obs['from']='GSE156625'
SOLD027.obs['from']='GSE156625'
SOLD028.obs['from']='GSE156625'
SOLD029.obs['from']='GSE156625'
SOLD030.obs['from']='GSE156625'
SOLD031.obs['from']='GSE156625'
SOLD032.obs['from']='GSE156625'
SOLD033.obs['from']='GSE156625'
SOLD034.obs['from']='GSE156625'
SOLD035.obs['from']='GSE156625'
SOLD036.obs['from']='GSE156625'
SOLD037.obs['from']='GSE156625'
SOLD038.obs['from']='GSE156625'
SOLD039.obs['from']='GSE156625'
SOLD040.obs['from']='GSE156625'
SOLD041.obs['from']='GSE156625'
SOLD042.obs['from']='GSE156625'
SOLD043.obs['from']='CRA002308'
SOLD044.obs['from']='CRA002308'
SOLD045.obs['from']='CRA002308'
SOLD046.obs['from']='CRA002308'
SOLD047.obs['from']='CRA002308'
SOLD048.obs['from']='CRA002308'
SOLD049.obs['from']='CRA002308'
SOLD050.obs['from']='CRA002308'
SOLD051.obs['from']='CRA002308'
SOLD052.obs['from']='CRA002308'
SOLD053.obs['from']='CRA002308'
SOLD054.obs['from']='CRA002308'
SOLD055.obs['from']='SRP318499'
SOLD056.obs['from']='SRP318499'
SOLD057.obs['from']='SRP318499'
SOLD058.obs['from']='SRP318499'
SOLD059.obs['from']='SRP318499'
SOLD060.obs['from']='SRP318499'
SOLD061.obs['from']='SRP318499'
SOLD062.obs['from']='GSE182159'
SOLD063.obs['from']='GSE182159'
SOLD064.obs['from']='GSE182159'
SOLD065.obs['from']='GSE182159'
SOLD066.obs['from']='GSE182159'
SOLD067.obs['from']='GSE182159'
SOLD068.obs['from']='GSE182159'
SOLD069.obs['from']='GSE182159'
SOLD070.obs['from']='GSE182159'
SOLD071.obs['from']='GSE182159'
SOLD072.obs['from']='GSE182159'
SOLD073.obs['from']='GSE182159'
SOLD074.obs['from']='GSE182159'
SOLD075.obs['from']='GSE182159'
SOLD076.obs['from']='GSE182159'
SOLD077.obs['from']='GSE182159'
SOLD078.obs['from']='GSE182159'
SOLD079.obs['from']='GSE182159'
SOLD080.obs['from']='GSE182159'
SOLD081.obs['from']='GSE182159'
SOLD082.obs['from']='GSE182159'
SOLD083.obs['from']='GSE182159'
SOLD084.obs['from']='GSE182159'
SOLD085.obs['from']='GSE252863'
SOLD086.obs['from']='GSE252863'
SOLD087.obs['from']='GSE252863'
SOLD088.obs['from']='GSE252863'
SOLD089.obs['from']='GSE252863'
SOLD090.obs['from']='GSE234241'
SOLD091.obs['from']='GSE234241'
SOLD092.obs['from']='GSE234241'
SOLD093.obs['from']='GSE234241'
SOLD094.obs['from']='GSE234241'
SOLD095.obs['from']='GSE234241'
SOLD096.obs['from']='GSE234241'
SOLD097.obs['from']='GSE234241'
SOLD098.obs['from']='GSE234241'
SOLD099.obs['from']='GSE234241'
SOLD100.obs['from']='GSE234241'
SOLD101.obs['from']='GSE186343'
SOLD102.obs['from']='GSE186343'
SOLD103.obs['from']='GSE186343'
SOLD104.obs['from']='GSE186343'
SOLD105.obs['from']='GSE186343'
SOLD106.obs['from']='GSE186343'




SOLD001.obs['Group'] = 'HCC'
SOLD002.obs['Group'] = 'HCC'
SOLD003.obs['Group'] = 'HCC'
SOLD004.obs['Group'] = 'HCC'
SOLD005.obs['Group'] = 'HCC'
SOLD006.obs['Group'] = 'HCC'
SOLD007.obs['Group'] = 'Normal'
SOLD008.obs['Group'] = 'Normal'
SOLD009.obs['Group'] = 'Normal'
SOLD010.obs['Group'] = 'Normal'
SOLD011.obs['Group'] = 'Normal'
SOLD012.obs['Group'] = 'Normal'
SOLD013.obs['Group'] = 'Normal'
SOLD014.obs['Group'] = 'Normal'
SOLD015.obs['Group'] = 'Normal'
SOLD016.obs['Group'] = 'HCC'
SOLD017.obs['Group'] = 'HCC'
SOLD018.obs['Group'] = 'HCC'
SOLD019.obs['Group'] = 'HCC'
SOLD020.obs['Group'] = 'HCC'
SOLD021.obs['Group'] = 'HCC'
SOLD022.obs['Group'] = 'HCC'
SOLD023.obs['Group'] = 'HCC'
SOLD024.obs['Group'] = 'HCC'
SOLD025.obs['Group'] = 'HCC'
SOLD026.obs['Group'] = 'HCC'
SOLD027.obs['Group'] = 'HCC'
SOLD028.obs['Group'] = 'HCC'
SOLD029.obs['Group'] = 'HCC'
SOLD030.obs['Group'] = 'HCC'
SOLD031.obs['Group'] = 'HCC'
SOLD032.obs['Group'] = 'HCC'
SOLD033.obs['Group'] = 'HCC'
SOLD034.obs['Group'] = 'HCC'
SOLD035.obs['Group'] = 'HCC'
SOLD036.obs['Group'] = 'HCC'
SOLD037.obs['Group'] = 'HCC'
SOLD038.obs['Group'] = 'HCC'
SOLD039.obs['Group'] = 'HCC'
SOLD040.obs['Group'] = 'HCC'
SOLD041.obs['Group'] = 'HCC'
SOLD042.obs['Group'] = 'HCC'
SOLD043.obs['Group'] = 'Normal'
SOLD044.obs['Group'] = 'HCC'
SOLD045.obs['Group'] = 'Normal'
SOLD046.obs['Group'] = 'HCC'
SOLD047.obs['Group'] = 'Normal'
SOLD048.obs['Group'] = 'HCC'
SOLD049.obs['Group'] = 'Normal'
SOLD050.obs['Group'] = 'HCC'
SOLD051.obs['Group'] = 'Normal'
SOLD052.obs['Group'] = 'HCC'
SOLD053.obs['Group'] = 'Normal'
SOLD054.obs['Group'] = 'HCC'
SOLD055.obs['Group'] = 'HCC'
SOLD056.obs['Group'] = 'HCC'
SOLD057.obs['Group'] = 'HCC'
SOLD058.obs['Group'] = 'HCC'
SOLD059.obs['Group'] = 'HCC'
SOLD060.obs['Group'] = 'HCC'
SOLD061.obs['Group'] = 'HCC'
SOLD062.obs['Group'] = 'Normal'
SOLD063.obs['Group'] = 'Normal'
SOLD064.obs['Group'] = 'Normal'
SOLD065.obs['Group'] = 'Normal'
SOLD066.obs['Group'] = 'Normal'
SOLD067.obs['Group'] = 'Normal'
SOLD068.obs['Group'] = 'Hepatitis'
SOLD069.obs['Group'] = 'Hepatitis'
SOLD070.obs['Group'] = 'Hepatitis'
SOLD071.obs['Group'] = 'Hepatitis'
SOLD072.obs['Group'] = 'Hepatitis'
SOLD073.obs['Group'] = 'Hepatitis'
SOLD074.obs['Group'] = 'Hepatitis'
SOLD075.obs['Group'] = 'Hepatitis'
SOLD076.obs['Group'] = 'Hepatitis'
SOLD077.obs['Group'] = 'Hepatitis'
SOLD078.obs['Group'] = 'Hepatitis'
SOLD079.obs['Group'] = 'Hepatitis'
SOLD080.obs['Group'] = 'Hepatitis'
SOLD081.obs['Group'] = 'Hepatitis'
SOLD082.obs['Group'] = 'Hepatitis'
SOLD083.obs['Group'] = 'Hepatitis'
SOLD084.obs['Group'] = 'Hepatitis'
SOLD085.obs['Group'] = 'Hepatitis'
SOLD086.obs['Group'] = 'Hepatitis'
SOLD087.obs['Group'] = 'Hepatitis'
SOLD088.obs['Group'] = 'Hepatitis'
SOLD089.obs['Group'] = 'Hepatitis'
SOLD090.obs['Group'] = 'Hepatitis'
SOLD091.obs['Group'] = 'Hepatitis'
SOLD092.obs['Group'] = 'Hepatitis'
SOLD093.obs['Group'] = 'Hepatitis'
SOLD094.obs['Group'] = 'Hepatitis'
SOLD095.obs['Group'] = 'Hepatitis'
SOLD096.obs['Group'] = 'Hepatitis'
SOLD097.obs['Group'] = 'Hepatitis'
SOLD098.obs['Group'] = 'Normal'
SOLD099.obs['Group'] = 'Normal'
SOLD100.obs['Group'] = 'Normal'
SOLD101.obs['Group'] = 'Hepatitis'
SOLD102.obs['Group'] = 'Hepatitis'
SOLD103.obs['Group'] = 'Hepatitis'
SOLD104.obs['Group'] = 'Hepatitis'
SOLD105.obs['Group'] = 'Hepatitis'
SOLD106.obs['Group'] = 'Hepatitis'


SOLD001.obs['Tissue'] = 'Tumor'
SOLD002.obs['Tissue'] = 'Tumor'
SOLD003.obs['Tissue'] = 'Tumor'
SOLD004.obs['Tissue'] = 'Tumor'
SOLD005.obs['Tissue'] = 'Tumor'
SOLD006.obs['Tissue'] = 'Tumor'
SOLD007.obs['Tissue'] = 'Normal'
SOLD008.obs['Tissue'] = 'Normal'
SOLD009.obs['Tissue'] = 'Normal'
SOLD010.obs['Tissue'] = 'Normal'
SOLD011.obs['Tissue'] = 'Normal'
SOLD012.obs['Tissue'] = 'Normal'
SOLD013.obs['Tissue'] = 'Normal'
SOLD014.obs['Tissue'] = 'Normal'
SOLD015.obs['Tissue'] = 'Normal'
SOLD016.obs['Tissue'] = 'Tumor'
SOLD017.obs['Tissue'] = 'Tumor'
SOLD018.obs['Tissue'] = 'Tumor'
SOLD019.obs['Tissue'] = 'Tumor'
SOLD020.obs['Tissue'] = 'Tumor'
SOLD021.obs['Tissue'] = 'Tumor'
SOLD022.obs['Tissue'] = 'Tumor'
SOLD023.obs['Tissue'] = 'Tumor'
SOLD024.obs['Tissue'] = 'Tumor'
SOLD025.obs['Tissue'] = 'Tumor'
SOLD026.obs['Tissue'] = 'Tumor'
SOLD027.obs['Tissue'] = 'Tumor'
SOLD028.obs['Tissue'] = 'Tumor'
SOLD029.obs['Tissue'] = 'Tumor'
SOLD030.obs['Tissue'] = 'Tumor'
SOLD031.obs['Tissue'] = 'Tumor'
SOLD032.obs['Tissue'] = 'Tumor'
SOLD033.obs['Tissue'] = 'Tumor'
SOLD034.obs['Tissue'] = 'Tumor'
SOLD035.obs['Tissue'] = 'Tumor'
SOLD036.obs['Tissue'] = 'Tumor'
SOLD037.obs['Tissue'] = 'Tumor'
SOLD038.obs['Tissue'] = 'Tumor'
SOLD039.obs['Tissue'] = 'Tumor'
SOLD040.obs['Tissue'] = 'Tumor'
SOLD041.obs['Tissue'] = 'Tumor'
SOLD042.obs['Tissue'] = 'Tumor'
SOLD043.obs['Tissue'] = 'Normal'
SOLD044.obs['Tissue'] = 'Tumor'
SOLD045.obs['Tissue'] = 'Normal'
SOLD046.obs['Tissue'] = 'Tumor'
SOLD047.obs['Tissue'] = 'Normal'
SOLD048.obs['Tissue'] = 'Tumor'
SOLD049.obs['Tissue'] = 'Normal'
SOLD050.obs['Tissue'] = 'Tumor'
SOLD051.obs['Tissue'] = 'Normal'
SOLD052.obs['Tissue'] = 'Tumor'
SOLD053.obs['Tissue'] = 'Normal'
SOLD054.obs['Tissue'] = 'Tumor'
SOLD055.obs['Tissue'] = 'Tumor'
SOLD056.obs['Tissue'] = 'Tumor'
SOLD057.obs['Tissue'] = 'Tumor'
SOLD058.obs['Tissue'] = 'Tumor'
SOLD059.obs['Tissue'] = 'Tumor'
SOLD060.obs['Tissue'] = 'Tumor'
SOLD061.obs['Tissue'] = 'Tumor'
SOLD062.obs['Tissue'] = 'Normal'
SOLD063.obs['Tissue'] = 'Normal'
SOLD064.obs['Tissue'] = 'Normal'
SOLD065.obs['Tissue'] = 'Normal'
SOLD066.obs['Tissue'] = 'Normal'
SOLD067.obs['Tissue'] = 'Normal'
SOLD068.obs['Tissue'] = 'Tumor'
SOLD069.obs['Tissue'] = 'Tumor'
SOLD070.obs['Tissue'] = 'Tumor'
SOLD071.obs['Tissue'] = 'Tumor'
SOLD072.obs['Tissue'] = 'Tumor'
SOLD073.obs['Tissue'] = 'Tumor'
SOLD074.obs['Tissue'] = 'Tumor'
SOLD075.obs['Tissue'] = 'Tumor'
SOLD076.obs['Tissue'] = 'Tumor'
SOLD077.obs['Tissue'] = 'Tumor'
SOLD078.obs['Tissue'] = 'Tumor'
SOLD079.obs['Tissue'] = 'Tumor'
SOLD080.obs['Tissue'] = 'Tumor'
SOLD081.obs['Tissue'] = 'Tumor'
SOLD082.obs['Tissue'] = 'Tumor'
SOLD083.obs['Tissue'] = 'Tumor'
SOLD084.obs['Tissue'] = 'Tumor'
SOLD085.obs['Tissue'] = 'Tumor'
SOLD086.obs['Tissue'] = 'Tumor'
SOLD087.obs['Tissue'] = 'Tumor'
SOLD088.obs['Tissue'] = 'Tumor'
SOLD089.obs['Tissue'] = 'Tumor'
SOLD090.obs['Tissue'] = 'Tumor'
SOLD091.obs['Tissue'] = 'Tumor'
SOLD092.obs['Tissue'] = 'Tumor'
SOLD093.obs['Tissue'] = 'Tumor'
SOLD094.obs['Tissue'] = 'Tumor'
SOLD095.obs['Tissue'] = 'Tumor'
SOLD096.obs['Tissue'] = 'Tumor'
SOLD097.obs['Tissue'] = 'Tumor'
SOLD098.obs['Tissue'] = 'Normal'
SOLD099.obs['Tissue'] = 'Normal'
SOLD100.obs['Tissue'] = 'Normal'
SOLD101.obs['Tissue'] = 'Tumor'
SOLD102.obs['Tissue'] = 'Tumor'
SOLD103.obs['Tissue'] = 'Tumor'
SOLD104.obs['Tissue'] = 'Tumor'
SOLD105.obs['Tissue'] = 'Tumor'
SOLD106.obs['Tissue'] = 'Tumor'


SOLD001.obs['SampleID']='SOLD001'
SOLD002.obs['SampleID']='SOLD002'
SOLD003.obs['SampleID']='SOLD003'
SOLD004.obs['SampleID']='SOLD004'
SOLD005.obs['SampleID']='SOLD005'
SOLD006.obs['SampleID']='SOLD006'
SOLD007.obs['SampleID']='SOLD007'
SOLD008.obs['SampleID']='SOLD008'
SOLD009.obs['SampleID']='SOLD009'
SOLD010.obs['SampleID']='SOLD010'
SOLD011.obs['SampleID']='SOLD011'
SOLD012.obs['SampleID']='SOLD012'
SOLD013.obs['SampleID']='SOLD013'
SOLD014.obs['SampleID']='SOLD014'
SOLD015.obs['SampleID']='SOLD015'
SOLD016.obs['SampleID']='SOLD016'
SOLD017.obs['SampleID']='SOLD017'
SOLD018.obs['SampleID']='SOLD018'
SOLD019.obs['SampleID']='SOLD019'
SOLD020.obs['SampleID']='SOLD020'
SOLD021.obs['SampleID']='SOLD021'
SOLD022.obs['SampleID']='SOLD022'
SOLD023.obs['SampleID']='SOLD023'
SOLD024.obs['SampleID']='SOLD024'
SOLD025.obs['SampleID']='SOLD025'
SOLD026.obs['SampleID']='SOLD026'
SOLD027.obs['SampleID']='SOLD027'
SOLD028.obs['SampleID']='SOLD028'
SOLD029.obs['SampleID']='SOLD029'
SOLD030.obs['SampleID']='SOLD030'
SOLD031.obs['SampleID']='SOLD031'
SOLD032.obs['SampleID']='SOLD032'
SOLD033.obs['SampleID']='SOLD033'
SOLD034.obs['SampleID']='SOLD034'
SOLD035.obs['SampleID']='SOLD035'
SOLD036.obs['SampleID']='SOLD036'
SOLD037.obs['SampleID']='SOLD037'
SOLD038.obs['SampleID']='SOLD038'
SOLD039.obs['SampleID']='SOLD039'
SOLD040.obs['SampleID']='SOLD040'
SOLD041.obs['SampleID']='SOLD041'
SOLD042.obs['SampleID']='SOLD042'
SOLD043.obs['SampleID']='SOLD043'
SOLD044.obs['SampleID']='SOLD044'
SOLD045.obs['SampleID']='SOLD045'
SOLD046.obs['SampleID']='SOLD046'
SOLD047.obs['SampleID']='SOLD047'
SOLD048.obs['SampleID']='SOLD048'
SOLD049.obs['SampleID']='SOLD049'
SOLD050.obs['SampleID']='SOLD050'
SOLD051.obs['SampleID']='SOLD051'
SOLD052.obs['SampleID']='SOLD052'
SOLD053.obs['SampleID']='SOLD053'
SOLD054.obs['SampleID']='SOLD054'
SOLD055.obs['SampleID']='SOLD055'
SOLD056.obs['SampleID']='SOLD056'
SOLD057.obs['SampleID']='SOLD057'
SOLD058.obs['SampleID']='SOLD058'
SOLD059.obs['SampleID']='SOLD059'
SOLD060.obs['SampleID']='SOLD060'
SOLD061.obs['SampleID']='SOLD061'
SOLD062.obs['SampleID']='SOLD062'
SOLD063.obs['SampleID']='SOLD063'
SOLD064.obs['SampleID']='SOLD064'
SOLD065.obs['SampleID']='SOLD065'
SOLD066.obs['SampleID']='SOLD066'
SOLD067.obs['SampleID']='SOLD067'
SOLD068.obs['SampleID']='SOLD068'
SOLD069.obs['SampleID']='SOLD069'
SOLD070.obs['SampleID']='SOLD070'
SOLD071.obs['SampleID']='SOLD071'
SOLD072.obs['SampleID']='SOLD072'
SOLD073.obs['SampleID']='SOLD073'
SOLD074.obs['SampleID']='SOLD074'
SOLD075.obs['SampleID']='SOLD075'
SOLD076.obs['SampleID']='SOLD076'
SOLD077.obs['SampleID']='SOLD077'
SOLD078.obs['SampleID']='SOLD078'
SOLD079.obs['SampleID']='SOLD079'
SOLD080.obs['SampleID']='SOLD080'
SOLD081.obs['SampleID']='SOLD081'
SOLD082.obs['SampleID']='SOLD082'
SOLD083.obs['SampleID']='SOLD083'
SOLD084.obs['SampleID']='SOLD084'
SOLD085.obs['SampleID']='SOLD085'
SOLD086.obs['SampleID']='SOLD086'
SOLD087.obs['SampleID']='SOLD087'
SOLD088.obs['SampleID']='SOLD088'
SOLD089.obs['SampleID']='SOLD089'
SOLD090.obs['SampleID']='SOLD090'
SOLD091.obs['SampleID']='SOLD091'
SOLD092.obs['SampleID']='SOLD092'
SOLD093.obs['SampleID']='SOLD093'
SOLD094.obs['SampleID']='SOLD094'
SOLD095.obs['SampleID']='SOLD095'
SOLD096.obs['SampleID']='SOLD096'
SOLD097.obs['SampleID']='SOLD097'
SOLD098.obs['SampleID']='SOLD098'
SOLD099.obs['SampleID']='SOLD099'
SOLD100.obs['SampleID']='SOLD100'
SOLD101.obs['SampleID']='SOLD101'
SOLD102.obs['SampleID']='SOLD102'
SOLD103.obs['SampleID']='SOLD103'
SOLD104.obs['SampleID']='SOLD104'
SOLD105.obs['SampleID']='SOLD105'
SOLD106.obs['SampleID']='SOLD106'


all_data_list = [SOLD001,SOLD002,SOLD003,SOLD004,SOLD005,SOLD006,SOLD007,SOLD008,SOLD009,SOLD010,SOLD011,SOLD012,SOLD013,SOLD014,SOLD015,SOLD016,SOLD017,SOLD018,SOLD019,SOLD020,SOLD021,SOLD022,SOLD023,SOLD024,SOLD025,SOLD026,SOLD027,SOLD028,SOLD029,SOLD030,SOLD031,SOLD032,SOLD033,SOLD034,SOLD035,SOLD036,SOLD037,SOLD038,SOLD039,SOLD040,SOLD041,SOLD042,SOLD043,SOLD044,SOLD045,SOLD046,SOLD047,SOLD048,SOLD049,SOLD050,SOLD051,SOLD052,SOLD053,SOLD054,SOLD055,SOLD056,SOLD057,SOLD058,SOLD059,SOLD060,SOLD061,SOLD062,SOLD063,SOLD064,SOLD065,SOLD066,SOLD067,SOLD068,SOLD069,SOLD070,SOLD071,SOLD072,SOLD073,SOLD074,SOLD075,SOLD076,SOLD077,SOLD078,SOLD079,SOLD080,SOLD081,SOLD082,SOLD083,SOLD084,SOLD085,SOLD086,SOLD087,SOLD088,SOLD089,SOLD090,SOLD091,SOLD092,SOLD093,SOLD094,SOLD095,SOLD096,SOLD097,SOLD098,SOLD099,SOLD100,SOLD101,SOLD102,SOLD103,SOLD104,SOLD105,SOLD106]
all_data = ad.concat(all_data_list, join="outer")


all_data.var_names_make_unique()
all_data.obs.drop(
    columns=['louvain','patient_id', 'patientno', 'patient_tumorsection', 'ViralvsNonViral',
             'PNC', 'PIC', 'NormalvsTumor', 'batch', 'batch_time',
             'sampleID'], 
    inplace=True)


new_var_names = [gene_name for gene_name in all_data.var_names if not gene_name.startswith('ENSG')]
all_data = all_data[:, new_var_names]
all_data

all_data.var['mt'] = all_data.var_names.str.startswith('MT-')
all_data.var['ribo'] = all_data.var_names.str.startswith(("RPS","RPL"))
sc.pp.calculate_qc_metrics(all_data, qc_vars=['mt','ribo'], percent_top=None, log1p=False, inplace=True)


sc.pl.violin(all_data, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter = 0.4, groupby = 'from')


all_data = all_data[all_data.obs.n_genes_by_counts < 6000, :]
all_data = all_data[all_data.obs.pct_counts_mt < 20, :]
sc.pp.filter_genes(all_data, min_cells = 3)
sc.pp.filter_cells(all_data, min_genes=200)

mito_genes = all_data.var_names.str.startswith('MT-')
hb_genes = all_data.var_names.str.contains('^HB[^(P)]')
ribo_genes = all_data.var_names.str.startswith(("RPS","RPL"))

remove = np.add(mito_genes, ribo_genes)
remove = np.add(remove, hb_genes)
keep = np.invert(remove)
all_data = all_data[:,keep]
all_data



