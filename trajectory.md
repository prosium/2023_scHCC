

## Preprocessing
```
library(devtools)
devtools::install_github("cole-trapnell-lab/monocle-release")

library(Seurat)
library(monocle)
library("SCOPfunctions")

load('./immune.combined_newAnno.singleR.RData')

data = as(utils_big_as.matrix(immune.combined.anno@assays$RNA@data, n_slices_init = 500, verbose = T), 'sparseMatrix')

pd = new("AnnotatedDataFrame", data = immune.combined.anno@meta.data)

newimport <- function(otherCDS, import_all = FALSE) {
  if(class(otherCDS)[1] == 'Seurat') {
    requireNamespace("Seurat")
    data <- otherCDS@assays$RNA@counts

    if(class(data) == "data.frame") {
      data <- as(as.matrix(data), "sparseMatrix")
    }

    pd <- tryCatch( {
      pd <- new("AnnotatedDataFrame", data = otherCDS@meta.data)
      pd
    }, 
    #warning = function(w) { },
    error = function(e) { 
      pData <- data.frame(cell_id = colnames(data), row.names = colnames(data))
      pd <- new("AnnotatedDataFrame", data = pData)

      message("This Seurat object doesn't provide any meta data");
      pd
    })

    # remove filtered cells from Seurat
    if(length(setdiff(colnames(data), rownames(pd))) > 0) {
      data <- data[, rownames(pd)]  
    }

    fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
    fd <- new("AnnotatedDataFrame", data = fData)
    lowerDetectionLimit <- 0

    if(all(data == floor(data))) {
      expressionFamily <- negbinomial.size()
    } else if(any(data < 0)){
      expressionFamily <- uninormal()
    } else {
      expressionFamily <- tobit()
    }

    valid_data <- data[, row.names(pd)]

    monocle_cds <- newCellDataSet(data,
                                  phenoData = pd, 
                                  featureData = fd,
                                  lowerDetectionLimit=lowerDetectionLimit,
                                  expressionFamily=expressionFamily)

    if(import_all) {
      if("Monocle" %in% names(otherCDS@misc)) {
        otherCDS@misc$Monocle@auxClusteringData$seurat <- NULL
        otherCDS@misc$Monocle@auxClusteringData$scran <- NULL

        monocle_cds <- otherCDS@misc$Monocle
        mist_list <- otherCDS

      } else {
        # mist_list <- list(ident = ident) 
        mist_list <- otherCDS
      }
    } else {
      mist_list <- list()
    }

    if(1==1) {
      var.genes <- setOrderingFilter(monocle_cds, otherCDS@assays$RNA@var.features)

    }
    monocle_cds@auxClusteringData$seurat <- mist_list

  } else if (class(otherCDS)[1] == 'SCESet') {
    requireNamespace("scater")

    message('Converting the exprs data in log scale back to original scale ...')    
    data <- 2^otherCDS@assayData$exprs - otherCDS@logExprsOffset

    fd <- otherCDS@featureData
    pd <- otherCDS@phenoData
    experimentData = otherCDS@experimentData
    if("is.expr" %in% slotNames(otherCDS))
      lowerDetectionLimit <- otherCDS@is.expr
    else 
      lowerDetectionLimit <- 1

    if(all(data == floor(data))) {
      expressionFamily <- negbinomial.size()
    } else if(any(data < 0)){
      expressionFamily <- uninormal()
    } else {
      expressionFamily <- tobit()
    }

    if(import_all) {
      # mist_list <- list(iotherCDS@sc3,
      #                   otherCDS@reducedDimension)
      mist_list <- otherCDS 

    } else {
      mist_list <- list()
    }

    monocle_cds <- newCellDataSet(data,
                                  phenoData = pd, 
                                  featureData = fd,
                                  lowerDetectionLimit=lowerDetectionLimit,
                                  expressionFamily=expressionFamily)
    # monocle_cds@auxClusteringData$sc3 <- otherCDS@sc3
    # monocle_cds@auxOrderingData$scran <- mist_list

    monocle_cds@auxOrderingData$scran <- mist_list

  } else {
    stop('the object type you want to export to is not supported yet')
  }

  return(monocle_cds)
}



save(Cholangiocytes, file = 'Cholangiocytes.Rdata')


library(Seurat)
library(monocle)
library("SCOPfunctions")

load('./Cholangiocytes.Rdata')


data = as(utils_big_as.matrix(Cholangiocytes@assays$RNA@data, n_slices_init = 100, verbose = T), 'sparseMatrix')

pd = new("AnnotatedDataFrame", data =Cholangiocytes@meta.data)


if(length(setdiff(colnames(data), rownames(pd))) > 0) {
    data <- data[, rownames(pd)]  
    }

fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new("AnnotatedDataFrame", data = fData)

valid_data <- data[, row.names(pd)]


monocle_cds <- newCellDataSet(valid_data,
                                  phenoData = pd, 
                                  featureData = fd,
                                  lowerDetectionLimit=0.5,
                                  expressionFamily=negbinomial.size())


monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds = estimateDispersions(monocle_cds)


monocle_cds = detectGenes(monocle_cds, min_expr = 0.1)
print(head(fData(monocle_cds)))
```


# Heatmap

```
setwd("/data/prosium/04.scRNA_HCC/20220517/")
BEAM_res <- readRDS("Epithelial_BEAM.res.rds")
monocle_cds <- readRDS("Epithelial_monocle.cds.rds")


load("../20220425/immune.combined_newAnno.singleR.RData") # annotated


target = "Cholangiocytes"
subset_obj <- subset(x=immune.combined.anno, idents = target)


subset_obj <- RunUMAP(subset_obj, dims = 1:10)
DimPlot(subset_obj, reduction = "umap",
             label = F, pt.size = 0.5, group.by = "celltype")+ NoLegend()

genes <- c("SQSTM1","SOD2","CCL2", "IL18", 
           "G0S2","DUSP5","ETS2","VEGFA","ATP2B1",
           "CCND1",
           "CXCL2","FOS","JUN","DUSP1","EGR1","IL6ST","HES1","CCNL1","ATF3","IRS2","FOSB","IER2")
genes2 <- c("CEBPB","NFIC")
FeaturePlot(subset_obj, features = genes2, ncol = 5)


VlnPlot(subset_obj, features = c("SOD2","SQSTM1","CCND1"), 
        slot = "counts", log = TRUE, group.by = "celltype")

plot_cell_trajectory(monocle_cds, color_by = "Pseudotime")

plot_genes_in_pseudotime(monocle_cds, color_by = "celltype")

plot_pseudotime_heatmap(HSMM_myo[sig_gene_names,],
                        num_clusters = 3,
                        cores = 1,
                        show_rownames = T)

```
