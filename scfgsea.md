## fgsea code
```
library(patchwork)
library(dplyr)
library(fgsea)
llibrary(pheatmap)

er=read.table(file="Normal_HBV_B_cell.DEG.txt",sep="\t",header=T)
er$fcsign <- sign(er$avg_log2FC)
er$logP=-log10(er$p_val)
er$metric= er$logP/er$fcsign
final<-er[,c("Gene","metric")]
pathways.disease <- gmtPathways("./h.all.v7.5.1.symbols.gmt")
write.table(er,file="./Normal_HBV_B_cell_GSEA.rnk",quote=F,sep="\t",row.names=F)
ranks <- deframe(final)
fgseaRes <- fgsea(pathways=pathways.disease, stats=ranks, nperm=1000)
fwrite(fgseaRes, file="./Normal_HBV_B_cell_fgseaRes.tsv", sep="\t", sep2=c("", " ", ""))


pval<-read.table(file="clipboard",sep="\t",header=T,row.names=1)
nes<-read.table(file="clipboard", sep="\t",header=T,row.names=1)

clab=matrix(pval, nc=ncol(pval), nr=nrow(pval))
clab[pval>0.05]=""
clab[pval<0.05]="*"
clab[pval<0.01]="**"
 pheatmap(nes, border_color="white", angle_col="315", fontsize_col=10, display_numbers=clab,number_color="black",fontsize_number=12, cellheight=10, cluster_cols=FALSE)

```
