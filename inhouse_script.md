
## 1. Prepare BAMs
```
wget https://sra-download.ncbi.nlm.nih.gov/traces/sra78/SRZ/014424/SRR14424777/740_possorted_genome_bam.bam
wget https://sra-download.ncbi.nlm.nih.gov/traces/sra34/SRZ/014424/SRR14424778/725_possorted_genome_bam.bam
wget https://sra-download.ncbi.nlm.nih.gov/traces/sra50/SRZ/014424/SRR14424779/713_possorted_genome_bam.bam
wget https://sra-download.ncbi.nlm.nih.gov/traces/sra1/SRZ/014424/SRR14424780/119_possorted_genome_bam.bam
wget https://sra-download.ncbi.nlm.nih.gov/traces/sra1/SRZ/014424/SRR14424781/114_possorted_genome_bam.bam
wget https://sra-download.ncbi.nlm.nih.gov/traces/sra33/SRZ/014424/SRR14424782/106_possorted_genome_bam.bam
wget https://sra-download.ncbi.nlm.nih.gov/traces/sra45/SRZ/014424/SRR14424783/104_possorted_genome_bam.bam
wget https://sra-download.ncbi.nlm.nih.gov/traces/sra60/SRZ/014424/SRR14424784/095_possorted_genome_bam.bam
```

## 2. BAMs to FASTQs
```
for file in $(ls *possorted*.bam);
do
    name=`echo -e ${file} | cut -f1 -d'.'`
    cellranger bamtofastq \
        --nthreads=8 \
        /mg5-90t/prj/2201_BSJ_HCC_scRNA/00.RAWDATA/${name}.bam \
        /mg5-90t/prj/2201_BSJ_HCC_scRNA/00.RAWDATA/${name}/
done
```


## 3. Count from FASTQs
```
for file in $(ls *possorted*/*/*I1_001.fastq.gz);
do
    name=`echo -e ${file} | cut -f1 -d'/'`
    path=`echo -e ${file} | cut -f1,2 -d'/'`
    echo -e ${name}
    echo -e ${PWD}/${path}
   
    cellranger count --id=${name}_count \
        --transcriptome=/mg5-90t/prj/2201_BSJ_HCC_scRNA/REF/refdata-gex-GRCh38-2020-A \
        --fastqs=${PWD}/${path} \
        --sample=${name} \
        --expect-cells=1000 \
        --localcores=16 \
        --localmem=100

done
```

## 4. Aggregation using Cellranger
```
touch library_info.csv
echo -e "sample_id,molecule_h5" >> library_info.csv
for file in $(ls /mg5-90t/prj/2201_BSJ_HCC_scRNA/01.HCC/00.RAWDATA/HCC*T/*_molecule_info.h5);
do
    name=`echo -e ${file} | cut -f7 -d'/'`
    echo -e "${name},${file}" >> library_info.csv
```

```
cellranger aggr \
    --id=HCC \
    --csv=library_info.csv \
    --normalize=mapped \
    --localcores=16 \
    --localmem=20
```

## 5. Quality control
```

library(dplyr)
library(Seurat)
library(patchwork)

# Load the PBMC dataset
aggr.data <- Read10X(data.dir = "./")
# Initialize the Seurat object with the raw (non-normalized data).
aggr <- CreateSeuratObject(counts = aggr.data, project = "HCC", min.cells = 3, min.features = 200)
> aggr
An object of class Seurat
22439 features across 30783 samples within 1 assay
Active assay: RNA (22439 features, 0 variable features)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
aggr[["percent.mt"]] <- PercentageFeatureSet(aggr, pattern = "^MT-")

tiff(filename = "01.HCC.QC.Vlnplot.tiff")
```


```
# Visualize QC metrics as a violin plot
VlnPlot(aggr, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
```
![image (1)](https://github.com/prosium/2023_scHCC/assets/94524627/1dd3eeaa-b520-4b1a-9ba7-8e5aadb5dfbf)


```
tiff(filename = "01.HCC.QC.FeatureScatter.tiff", width = 600, height = 600)
plot1 <- FeatureScatter(aggr, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(aggr, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()
```
![image (2)](https://github.com/prosium/2023_scHCC/assets/94524627/1342c264-aac4-4c01-a90b-16c58e38d95a)


## 6. Filtering

```
aggr <- subset(aggr, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

> aggr
An object of class Seurat
22439 features across 15051 samples within 1 assay
Active assay: RNA (22439 features, 0 variable features)
```

## 7. Normalization

```
aggr <- NormalizeData(aggr, normalization.method = "LogNormalize", scale.factor = 10000)
```


## 8. Identification of highly varaible features

```
aggr <- FindVariableFeatures(aggr, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(aggr), 10)

tiff(filename = "./02.HCC.variable_feature.tiff", width = 600, height = 600)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(aggr)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
dev.off()
```
![image (3)](https://github.com/prosium/2023_scHCC/assets/94524627/048b1ceb-822d-4eea-8c01-232a1787f81c)


## 9. Scaling
```
all.genes <- rownames(aggr)
aggr <- ScaleData(aggr, features = all.genes)
```

## 10. Linear dimension reduction
```
aggr <- RunPCA(aggr, features = VariableFeatures(object = aggr))

> print(aggr[["pca"]], dims = 1:5, nfeatures = 5)
PC_ 1
Positive:  SPARC, IGFBP7, GNG11, COL4A1, FLT1
Negative:  APOC1, APOA2, APOH, APOC3, RBP4
PC_ 2
Positive:  ACTA2, RGS5, TAGLN, CALD1, MYL9
Negative:  C1QB, C1QA, C1QC, HLA-DRA, HLA-DRB1
PC_ 3
Positive:  APOC2, APOH, SERPINA1, RBP4, AMBP
Negative:  KLRB1, NKG7, LTB, CD69, ACTA2
PC_ 4
Positive:  IL32, KLRB1, RGCC, NKG7, LTB
Negative:  TAGLN, ACTA2, RGS5, MYL9, TIMP1
PC_ 5
Positive:  FXYD2, SNORC, PTGDS, TTYH1, ANXA13
Negative:  CYP2E1, HPD, CFHR2, FGB, GSTA1
```


## 11. Determine the dimensionality
```
tiff('./03.HCC.JackStrawPlot.tiff')
JackStrawPlot(aggr, dims = 1:20)
dev.off()
```
![image (4)](https://github.com/prosium/2023_scHCC/assets/94524627/a6d5135b-6bed-4e16-b369-e60a67596412)

```
tiff('./03.HCC.ElbowPlot.tiff')
ElbowPlot(aggr)
dev.off()
```
![image (5)](https://github.com/prosium/2023_scHCC/assets/94524627/8e903edc-8901-40f1-a63d-d5e576102a19)


## 12. Cluster cells
```
aggr <- FindNeighbors(aggr, dims = 1:15)
aggr <- FindClusters(aggr, resolution = 0.5)

> head(Idents(aggr), 5)
AAACCCACAATAAGGT-1 AAACCCACACGTGAGA-1 AAACCCACAGCTAACT-1 AAACCCACATCAACCA-1
                 9                  9                  9                 14
AAACCCAGTTACGTAC-1
                 9
Levels: 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17

aggr <- RunUMAP(aggr, dims = 1:15)

pdf('./04.HCC.cluster.UMAP.pdf')
DimPlot(aggr, reduction = "umap", label = TRUE)
dev.off()
```
![image (6)](https://github.com/prosium/2023_scHCC/assets/94524627/b6acfbbc-981c-4a95-b39d-2cd103202bca)


## 13. Cluster cells
```
library(limma)
cluster = strsplit2(rownames(aggr@meta.data), '\\-')[,2]
cluster = as.matrix(factor(cluster, label = c('HCC10T', 'HCC11T', 'HCC12T', 'HCC2T', 'HCC3T', 'HCC5T', 'HCC8T', 'HCC9T')))
rownames(cluster) = rownames(aggr@meta.data)

# add some more meta data
aggr <- AddMetaData(object = aggr,
                    metadata = cluster,
                    col.name = "cluster")
```
![image (7)](https://github.com/prosium/2023_scHCC/assets/94524627/b5de9123-54f4-4fa4-8e02-82946da1860c)

## 14. SingleR annotation
```

library(SingleR)
library(dplyr)
library(patchwork)
library(MASS)
library(celldex)
library(SingleCellExperiment)
hpca.se <- HumanPrimaryCellAtlasData()

results <- SingleR(test = as.SingleCellExperiment(aggr),
                   ref = list(hpca.se),
                   labels = list(hpca.se$label.main))

# Copy over the labels and pruned.labels
aggr$SingleR.pruned.calls <- results$pruned.labels
aggr$SingleR.calls <- results$labels

write.table(aggr@meta.data, file = "./HCC_aggr.results_anno.txt", sep = '\t')
```














