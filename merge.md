## merge code
```
library(patchwork)
library(Seurat)
GSE151530<-Read10X(data.dir="D:/sjbaek/test/GSE151530")
GSE151530<-CreateSeuratObject(counts=GSE151530, project="GSE151530")
metasample=read.table("D:/sjbaek/test/GSE151530_Clin.txt",sep="\t",header=T)
head(metasample)
    ID Sample               cell    type diagnosis stage etiology treatment
1 S001    H08 AAACCTGAGCTCTCGG-1 T cells  American   HCC       IV       HCV
2 S001    H08 AAACCTGCAAGAAAGG-1 T cells  American   HCC       IV       HCV
3 S001    H08 AAACCTGCAAGTACCT-1 T cells  American   HCC       IV       HCV
4 S001    H08 AAACCTGCAAGTAGTA-1 B cells  American   HCC       IV       HCV
5 S001    H08 AAACCTGGTCGGATCC-1    TAMs  American   HCC       IV       HCV
6 S001    H08 AAACCTGGTTACGGAG-1    TAMs  American   HCC       IV       HCV
> metasample$virus=metasample$treatment

GSE151530@meta.data$Cell=rownames(GSE151530@meta.data)
head(GSE151530@meta.data)
                   orig.ident nCount_RNA nFeature_RNA               Cell
AAACCTGAGCTCTCGG-1  GSE151530       1688          646 AAACCTGAGCTCTCGG-1
AAACCTGCAAGAAAGG-1  GSE151530       4352          645 AAACCTGCAAGAAAGG-1
AAACCTGCAAGTACCT-1  GSE151530       1794          674 AAACCTGCAAGTACCT-1
AAACCTGCAAGTAGTA-1  GSE151530       2953          848 AAACCTGCAAGTAGTA-1
AAACCTGGTCGGATCC-1  GSE151530       1902          850 AAACCTGGTCGGATCC-1
AAACCTGGTTACGGAG-1  GSE151530       4061         1377 AAACCTGGTTACGGAG-1

library(rlang)
meta_data<-dplyr::left_join(GSE151530@meta.data, metasample, by="Cell")

> head(metasample)
    ID Sample               cell    type diagnosis stage etiology treatment virus
1 S001    H08 AAACCTGAGCTCTCGG-1 T cells  American   HCC       IV       HCV   HCV
2 S001    H08 AAACCTGCAAGAAAGG-1 T cells  American   HCC       IV       HCV   HCV
3 S001    H08 AAACCTGCAAGTACCT-1 T cells  American   HCC       IV       HCV   HCV
4 S001    H08 AAACCTGCAAGTAGTA-1 B cells  American   HCC       IV       HCV   HCV
5 S001    H08 AAACCTGGTCGGATCC-1    TAMs  American   HCC       IV       HCV   HCV
6 S001    H08 AAACCTGGTTACGGAG-1    TAMs  American   HCC       IV       HCV   HCV
GSE151530@meta.data$cell=rownames(GSE151530@meta.data)
meta_data<-dplyr::left_join(GSE151530@meta.data, metasample, by="cell")
> head(meta_data)
  orig.ident nCount_RNA nFeature_RNA               Cell               cell   ID Sample    type diagnosis stage etiology treatment virus
1  GSE151530       1688          646 AAACCTGAGCTCTCGG-1 AAACCTGAGCTCTCGG-1 S001    H08 T cells  American   HCC       IV       HCV   HCV
2  GSE151530       4352          645 AAACCTGCAAGAAAGG-1 AAACCTGCAAGAAAGG-1 S001    H08 T cells  American   HCC       IV       HCV   HCV
3  GSE151530       1794          674 AAACCTGCAAGTACCT-1 AAACCTGCAAGTACCT-1 S001    H08 T cells  American   HCC       IV       HCV   HCV
4  GSE151530       2953          848 AAACCTGCAAGTAGTA-1 AAACCTGCAAGTAGTA-1 S001    H08 B cells  American   HCC       IV       HCV   HCV
5  GSE151530       1902          850 AAACCTGGTCGGATCC-1 AAACCTGGTCGGATCC-1 S001    H08    TAMs  American   HCC       IV       HCV   HCV
6  GSE151530       4061         1377 AAACCTGGTTACGGAG-1 AAACCTGGTTACGGAG-1 S001    H08    TAMs  American   HCC       IV       HCV   HCV
 rownames(meta_data)=GSE151530@meta.data$cell

GSE151530@meta.data<-meta_data
 Idents(GSE151530)=GSE151530@meta.data$treatment
HBV=subset(x=GSE151530, idents="HBV")
HCV=subset(x=GSE151530, idents="HCV")
GSE151530_HBV<-HBV
GSE151530_HCV<-HCV

####SRP318499(HBV)
> SRR14424777<-Read10X(data.dir="D:/sjbaek/test/SRR14424777")
> SRR14424777<-CreateSeuratObject(counts=SRR14424777, project="SRR14424777")
> SRR14424778<-Read10X(data.dir="D:/sjbaek/test/SRR14424778")
> SRR14424778<-CreateSeuratObject(counts=SRR14424778, project="SRR14424778")
> SRR14424779<-Read10X(data.dir="D:/sjbaek/test/SRR14424779")
> SRR14424779<-CreateSeuratObject(counts=SRR14424779, project="SRR14424779")
> SRR14424781<-Read10X(data.dir="D:/sjbaek/test/SRR14424781")
> SRR14424781<-CreateSeuratObject(counts=SRR14424781, project="SRR14424781")
> SRR14424782<-Read10X(data.dir="D:/sjbaek/test/SRR14424782")
> SRR14424782<-CreateSeuratObject(counts=SRR14424782, project="SRR14424782")
> SRR14424783<-Read10X(data.dir="D:/sjbaek/test/SRR14424783")
> SRR14424783<-CreateSeuratObject(counts=SRR14424783, project="SRR14424783")
> SRR14424784<-Read10X(data.dir="D:/sjbaek/test/SRR14424784")
> SRR14424784<-CreateSeuratObject(counts=SRR14424784, project="SRR14424784")
> SRP318499_HBV<-merge(SRR14424777, y=list(SRR14424778, SRR14424779, SRR14424781, SRR14424782, SRR14424783, SRR14424784), add.cell.ids=c("SRR14424777","SRR14424778","SRR14424779","SRR14424781","SRR14424782","SRR14424783","SRR14424784"),project="HBV")

> head(SRP318499_HBV 7 )
                                orig.ident nCount_RNA nFeature_RNA
SRR14424777_AAACCTGAGAAACCAT-1 SRR14424777        427          280
SRR14424777_AAACCTGAGAAACCGC-1 SRR14424777        219          162
SRR14424777_AAACCTGAGAAACCTA-1 SRR14424777        202          144
SRR14424777_AAACCTGAGAAACGAG-1 SRR14424777          5            5
SRR14424777_AAACCTGAGAAACGCC-1 SRR14424777          2            2
SRR14424777_AAACCTGAGAAAGTGG-1 SRR14424777          5            4
SRR14424777_AAACCTGAGAACAACT-1 SRR14424777          3            3
SRR14424777_AAACCTGAGAACAATC-1 SRR14424777          3            3
SRR14424777_AAACCTGAGAACTCGG-1 SRR14424777          4            4
SRR14424777_AAACCTGAGAAGAAGC-1 SRR14424777          2            2

####################CRA002308 (HBV 6 ; HCV; 1)
> H01T<-Read10X(data.dir="D:/sjbaek/test/CRR115704")
> H01T<-CreateSeuratObject(counts=H01T, project="H01T")
> H02T<-Read10X(data.dir="D:/sjbaek/test/CRR115706")
> H02T<-CreateSeuratObject(counts=H02T, project="H02T")
> H03T<-Read10X(data.dir="D:/sjbaek/test/CRR115708")
> H03T<-CreateSeuratObject(counts=H03T, project="H03T")
> H04T<-Read10X(data.dir="D:/sjbaek/test/CRR115710")
> H04T<-CreateSeuratObject(counts=H04T, project="H04T")
> H05T<-Read10X(data.dir="D:/sjbaek/test/CRR115711")
> H05T<-CreateSeuratObject(counts=H05T, project="H05T")
> H06T<-Read10X(data.dir="D:/sjbaek/test/CRR115714")
> H06T<-CreateSeuratObject(counts=H06T, project="H06T")
> H07T<-Read10X(data.dir="D:/sjbaek/test/CRR115716")
> H07T<-CreateSeuratObject(counts=H07T, project="H07T")
> CRA002308_HBV<-merge(H01T, y=list(H02T, H03T, H04T, H06T, H07T), add.cell.ids=c("H01T","H02T","H03T","H04T","H06T","H07T"), project="HBV")
> CRA002308_HCV<-H05T
> save(CRA002308_HBV, file="D:/sjbaek/test/CRA002308_HBV.RData")
> save(CRA002308_HCV, file="D:/sjbaek/test/CRA002308_HCV.RData")
> HCC2T.data<-Read10X(data.dir="D:/sjbaek/test/HCC2T")
> HCC3T.data<-Read10X(data.dir="D:/sjbaek/test/HCC3T")
> HCC5T.data<-Read10X(data.dir="D:/sjbaek/test/HCC5T")
> HCC8T.data<-Read10X(data.dir="D:/sjbaek/test/HCC8T")
> HCC9T.data<-Read10X(data.dir="D:/sjbaek/test/HCC9T")
> HCC10T.data<-Read10X(data.dir="D:/sjbaek/test/HCC10T")
> HCC11T.data<-Read10X(data.dir="D:/sjbaek/test/HCC11T")
> HCC12T.data<-Read10X(data.dir="D:/sjbaek/test/HCC12T")
> HCC2T<-CreateSeuratObject(counts=HCC2T.data, project="HCC2T")
> HCC3T<-CreateSeuratObject(counts=HCC3T.data, project="HCC3T")
> HCC5T<-CreateSeuratObject(counts=HCC5T.data, project="HCC5T")
> HCC8T<-CreateSeuratObject(counts=HCC8T.data, project="HCC8T")
> HCC9T<-CreateSeuratObject(counts=HCC9T.data, project="HCC9T")
> HCC10T<-CreateSeuratObject(counts=HCC10T.data, project="HCC10T")
> HCC11T<-CreateSeuratObject(counts=HCC11T.data, project="HCC11T")
> HCC12T<-CreateSeuratObject(counts=HCC12T.data, project="HCC12T")
> home_HBV<-merge(HCC2T, y=list( HCC5T, HCC8T, HCC10T, HCC11T, HCC12T), add.cell.ids=c("HCC2T","HCC5T","HCC8T","HCC10T","HCC11T","HCC12T"), project="HBV")
> home_HCV<-merge(HCC3T, y=HCC9T, add.cell.ids=c("HCC3T","HCC9T"), project="HCV")
> save(home_HBV, file="D:/sjbaek/test/home_HBV.RData")
> save(home_HCV, file="D:/sjbaek/test/home_HCV.RData")
> head(home_HCV)
                         orig.ident nCount_RNA nFeature_RNA
HCC3T_AAACCTGAGACACGAC-1      HCC3T       3243         1270
HCC3T_AAACCTGAGAGATGAG-1      HCC3T        500          355
HCC3T_AAACCTGGTCAGATAA-1      HCC3T       2053          714
HCC3T_AAACGGGAGCTCCTTC-1      HCC3T       6158         1665
HCC3T_AAACGGGCAGGTGGAT-1      HCC3T      22811         3979
HCC3T_AAACGGGCATCGATTG-1      HCC3T        935          520
HCC3T_AAACGGGTCGCATGGC-1      HCC3T      13652         2962
HCC3T_AAAGATGAGATGCGAC-1      HCC3T       1272          670
HCC3T_AAAGATGAGGAGTCTG-1      HCC3T       1105          537
HCC3T_AAAGATGTCCACGACG-1      HCC3T      15628         3086

======================merge each HBV, HCV ================================
            
HBV_merge<-merge(home_HBV,y=list(SRP318499_HBV,CRA002308_HBV, GSE151530_HBV), add.cell.ids=c("homeHBV","SRP318499HBV","CRA002308HBV","GSE151530HBV"), project="HBV")
HCV_merge<-merge(home_HCV, y=list(CRA002308_HCV, GSE151530_HCV), add.cell.ids=c("homeHCV","CRA002308HCV","GSE151530HCV"),project="HCV")
HBV_HCV_merge<-merge(HBV_merge,y=list(HCV_merge), add.cell.ids=c("HBV","HCV"), project="HCC")

=========================================merge HBV and HCV (ALL)=================

> load("D:/KIOM/work/HCC_singlecell/public_merged/HCV_merge.RData")

> library(Seurat)

> load("D:/KIOM/work/HCC_singlecell/public_merged/HBV_merge.RData")

##merge

> HBV_HCV_merge<-merge(HBV_merge,y=list(HCV_merge), add.cell.ids=c("HBV","HCV"), project="HCC")



## batch correction

lifnb.list <- SplitObject(HBV_HCV_merge, split.by = "celltype")


ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = ifnb.list)
All.anchors<-FindIntegrationAnchors(object.list=ifnb.list, anchor.features=features)



All.combined<-IntegratedData(anchorset=immune.anchors)
```
