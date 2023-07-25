```
library(Seurat)

load('./immune.combined.RData')

library(SingleR)
library(scater)
library(dplyr)
library(patchwork)
library(MASS)

library(celldex)
library(SingleCellExperiment)
hpca.se <- HumanPrimaryCellAtlasData()

results <- SingleR(test = as.SingleCellExperiment(immune.combined),
                        ref = list(hpca.se),
                        labels = list(hpca.se$label.main))

# Copy over the labels and pruned.labels
immune.combined$SingleR.pruned.calls <- results$pruned.labels
immune.combined$SingleR.calls <- results$labels

write.table(immune.combined@meta.data, file = 'immune.combined.meta_data.txt', sep = '\t', quote = F)

### 각 클러스터별 annotation 통계 ###

#### cluster 0 ####
Cluster0 = subset(x = immune.combined, idents = '0')

> table(Cluster0$SingleR.calls)

     Hepatocytes          NK_cell Pre-B_cell_CD34-          T_cells
               4              407                2            16849

=> T_cells
#### cluster 1 ####
Cluster1 = subset(x = immune.combined, idents = '1')

> table(Cluster1$SingleR.calls)

          B_cell Epithelial_cells              GMP      Hepatocytes
               4                2                1            15248
        Monocyte          Neurons      Neutrophils          NK_cell
               1                1                2                4
Pro-B_cell_CD34+          T_cells
               3               14

=> Hepatocytes
#### cluster 2 ####
Cluster2 = subset(x = immune.combined, idents = '2')

> table(Cluster2$SingleR.calls)

Hepatocytes  HSC_-G-CSF     NK_cell     T_cells
          3           1         284       12830

=> T_cells
#### cluster 3 ####
Cluster3 = subset(x = immune.combined, idents = '3')

> table(Cluster3$SingleR.calls)

          B_cell     Erythroblast      Gametocytes              GMP
               6                4                4                1
     Hepatocytes        HSC_CD34+              MEP Pro-B_cell_CD34+
           10320                1                1                2
         T_cells
              20

=> Hepatocytes
#### cluster 4 ####
Cluster4 = subset(x = immune.combined, idents = '4')

> table(Cluster4$SingleR.calls)

NK_cell T_cells
   7628     194

=> NK_cell
#### cluster 5 ####
Cluster5 = subset(x = immune.combined, idents = '5')

> table(Cluster5$SingleR.calls)

         NK_cell Pre-B_cell_CD34-          T_cells
            7417                2              242

=> NK_cell
#### cluster 6 ####
Cluster6 = subset(x = immune.combined, idents = '6')

> table(Cluster6$SingleR.calls)

          B_cell     Erythroblast      Gametocytes      Hepatocytes
               3                1                3             7476
      Macrophage          NK_cell Pre-B_cell_CD34- Pro-B_cell_CD34+
               1                1                1                1
         T_cells
              19

=> Hepatocytes
#### cluster 7 ####
Cluster7 = subset(x = immune.combined, idents = '7')

> table(Cluster7$SingleR.calls)

              B_cell Embryonic_stem_cells     Epithelial_cells
                   4                    1                    1
         Gametocytes          Hepatocytes            HSC_CD34+
                   9                 6655                    1
           iPS_cells        Keratinocytes              NK_cell
                   1                    2                    3
    Pro-B_cell_CD34+              T_cells
                   1                   10

=> Hepatocytes
#### cluster 8 ####
Cluster8 = subset(x = immune.combined, idents = '8')

> table(Cluster8$SingleR.calls)

              DC      Hepatocytes       HSC_-G-CSF        HSC_CD34+
             132               67               84                1
      Macrophage         Monocyte      Neutrophils          NK_cell
             399             5786               27               11
Pre-B_cell_CD34-          T_cells
              30               33

=> Monocyte     
#### cluster 9 ####
Cluster9 = subset(x = immune.combined, idents = '9')

> table(Cluster9$SingleR.calls)

          B_cell              CMP      Hepatocytes       HSC_-G-CSF
               2                1               10                2
      Macrophage          NK_cell Pre-B_cell_CD34- Pro-B_cell_CD34+
               1               42                2                1
         T_cells
            6487

=> T_cells
#### cluster 10 ####

Cluster10 = subset(x = immune.combined, idents = '10')

> table(Cluster10$SingleR.calls)

NK_cell T_cells
   1671    4522

=> T_cells
#### cluster 11 ####

Cluster11 = subset(x = immune.combined, idents = '11')

> table(Cluster11$SingleR.calls)

NK_cell T_cells
    731    5378

=> T_cells
#### cluster 12 ####

Cluster12 = subset(x = immune.combined, idents = '12')

> table(Cluster12$SingleR.calls)

           Astrocyte               B_cell         Chondrocytes
                   3                    1                 1637
Embryonic_stem_cells    Endothelial_cells          Fibroblasts
                   7                    3                  537
         Hepatocytes           Macrophage             Monocyte
                 191                    1                    3
                 MSC Neuroepithelial_cell              Neurons
                  84                    7                    7
         Neutrophils              NK_cell          Osteoblasts
                   1                    3                   61
Smooth_muscle_cells              T_cells    Tissue_stem_cells
                2030                   16                 1402

=> Smooth_muscle_cells
#### cluster 13 ####

Cluster13 = subset(x = immune.combined, idents = '13')

> table(Cluster13$SingleR.calls)

           Astrocyte               B_cell         Chondrocytes
                  39                    3                  128
                  DC Embryonic_stem_cells    Endothelial_cells
                  11                    8                 3435
    Epithelial_cells          Fibroblasts                  GMP
                   4                  119                    1
         Hepatocytes           HSC_-G-CSF            HSC_CD34+
                1033                    3                    6
           iPS_cells           Macrophage             Monocyte
                   3                   19                    4
                 MSC Neuroepithelial_cell              Neurons
                  10                    4                   59
         Neutrophils              NK_cell          Osteoblasts
                   1                  104                   13
    Pre-B_cell_CD34-     Pro-B_cell_CD34+  Smooth_muscle_cells
                   2                    2                   86
             T_cells    Tissue_stem_cells
                 207                  190

=> Endothelial_cells
#### cluster 14 ####

Cluster14 = subset(x = immune.combined, idents = '14')

> table(Cluster14$SingleR.calls)

          B_cell              CMP               DC      Gametocytes
              25                4              668                2
             GMP      Hepatocytes       HSC_-G-CSF        HSC_CD34+
              40             1642               57                6
      Macrophage         Monocyte        Myelocyte      Neutrophils
            1235              730                1               19
         NK_cell        Platelets Pre-B_cell_CD34- Pro-B_cell_CD34+
             285                1              168               10
   Pro-Myelocyte          T_cells
               3              559

=> Hepatocytes
#### cluster 15 ####

Cluster15 = subset(x = immune.combined, idents = '15')

> table(Cluster15$SingleR.calls)

          B_cell              CMP               DC              GMP
            4777                2                1                1
     Hepatocytes       Macrophage         Monocyte          NK_cell
              59                1                1               28
Pro-B_cell_CD34+          T_cells
               6              273

=> B_cell             
#### cluster 16 ####

Cluster16 = subset(x = immune.combined, idents = '16')

> table(Cluster16$SingleR.calls)

        Astrocyte            B_cell               CMP                DC
                1                 7                 3                10
Endothelial_cells               GMP       Hepatocytes        HSC_-G-CSF
                2                 2               748                23
        HSC_CD34+        Macrophage               MEP          Monocyte
                5                 8                 1                24
      Neutrophils           NK_cell  Pre-B_cell_CD34-  Pro-B_cell_CD34+
               11              1878                24                 2
          T_cells
             2382

= > T_cells
#### cluster 17 ####

Cluster17 = subset(x = immune.combined, idents = '17')

> table(Cluster17$SingleR.calls)

         DC Hepatocytes   HSC_CD34+  Macrophage    Monocyte Neutrophils
        410         561           2        3737         290           6
    NK_cell     T_cells
         16          13

=> Macrophage   
#### cluster 18 ####

Cluster18 = subset(x = immune.combined, idents = '18')

> table(Cluster18$SingleR.calls)

Hepatocytes     NK_cell     T_cells
       1292          33        2806

=> T_cells
#### cluster 19 ####

Cluster19 = subset(x = immune.combined, idents = '19')

> table(Cluster19$SingleR.calls)

  Gametocytes   Hepatocytes Keratinocytes   Neutrophils
            1          4119             1             1

=> Hepatocytes
#### cluster 20 ####

Cluster20 = subset(x = immune.combined, idents = '20')

> table(Cluster20$SingleR.calls)

           Astrocyte               B_cell         Chondrocytes
                   6                    3                   87
                 CMP                   DC Embryonic_stem_cells
                   3                    1                   12
   Endothelial_cells     Epithelial_cells          Fibroblasts
                   9                    1                   13
                 GMP          Hepatocytes            iPS_cells
                   1                 2628                    1
          Macrophage             Monocyte                  MSC
                   2                    3                    4
Neuroepithelial_cell              Neurons              NK_cell
                  16                    3                   65
         Osteoblasts            Platelets     Pro-B_cell_CD34+
                   2                    1                    1
Smooth_muscle_cells              T_cells    Tissue_stem_cells
                  38                  417                   77

=> Hepatocytes           
#### cluster 21 ####

Cluster21 = subset(x = immune.combined, idents = '21')

> table(Cluster21$SingleR.calls)

Hepatocytes     NK_cell     T_cells
        756         390        2107

=> T_cells
#### cluster 22 ####

Cluster22 = subset(x = immune.combined, idents = '22')

> table(Cluster22$SingleR.calls)

     B_cell Gametocytes Hepatocytes    Monocyte     NK_cell     T_cells
          3           1        3194           1           3          40

=> Hepatocytes   
#### cluster 23 ####

Cluster23 = subset(x = immune.combined, idents = '23')

> table(Cluster23$SingleR.calls)

             B_cell    Epithelial_cells         Hepatocytes           Myelocyte
                  1                   2                3029                   1
            NK_cell Smooth_muscle_cells             T_cells
                  1                   1                  41

=> Hepatocytes
#### cluster 24 ####

Cluster24 = subset(x = immune.combined, idents = '24')

> table(Cluster24$SingleR.calls)

          B_cell              CMP               DC              GMP
            1368                2                1                2
     Hepatocytes        HSC_CD34+       Macrophage         Monocyte
             550                1                1                5
         NK_cell Pre-B_cell_CD34- Pro-B_cell_CD34+    Pro-Myelocyte
              25               10               42                1
         T_cells
             620

=> B_cell             
#### cluster 25 ####

Cluster25 = subset(x = immune.combined, idents = '25')

> table(Cluster25$SingleR.calls)

         Fibroblasts          Hepatocytes             Monocyte
                   1                 2598                    1
Neuroepithelial_cell              T_cells
                   1

=> Hepatocytes            
#### cluster 26 ####

Cluster26 = subset(x = immune.combined, idents = '26')

> table(Cluster26$SingleR.calls)

Hepatocytes
       2142

=> Hepatocytes
#### cluster 27 ####

Cluster27 = subset(x = immune.combined, idents = '27')

> table(Cluster27$SingleR.calls)

                  BM Embryonic_stem_cells     Epithelial_cells
                   1                    1                    2
         Hepatocytes            Myelocyte              Neurons
                2095                    1                    2
Smooth_muscle_cells              T_cells
                   1                   17

=> Hepatocytes           
#### cluster 28 ####

Cluster28 = subset(x = immune.combined, idents = '28')

> table(Cluster28$SingleR.calls)

              B_cell           BM & Prog.                  CMP
                  50                    2                    7
Embryonic_stem_cells          Hepatocytes           HSC_-G-CSF
                   2                  123                    1
          Macrophage                  MEP Neuroepithelial_cell
                   1                    1                    2
         Neutrophils              NK_cell     Pre-B_cell_CD34-
                   1                  265                    3
    Pro-B_cell_CD34+              T_cells
                  24                 1529

=> T_cells
#### cluster 29 ####

Cluster29 = subset(x = immune.combined, idents = '29')

> table(Cluster29$SingleR.calls)

Hepatocytes     NK_cell     T_cells
         41         794        1040

=> T_cells
#### cluster 30 ####

Cluster30 = subset(x = immune.combined, idents = '30')

> table(Cluster30$SingleR.calls)

Hepatocytes
       1872

=> Hepatocytes
#### cluster 31 ####

Cluster31 = subset(x = immune.combined, idents = '31')

> table(Cluster31$SingleR.calls)

Hepatocytes     NK_cell     T_cells
       1376         160         332

=> Hepatocytes    
#### cluster 32 ####

Cluster32 = subset(x = immune.combined, idents = '32')

> table(Cluster32$SingleR.calls)

Hepatocytes     NK_cell     T_cells
        128         821         722

=> NK_cell    
#### cluster 33 ####

Cluster33 = subset(x = immune.combined, idents = '33')

> table(Cluster33$SingleR.calls)

           Astrocyte           BM & Prog.         Chondrocytes
                   3                    2                   13
                 CMP                   DC Embryonic_stem_cells
                   2                    3                   12
   Endothelial_cells     Epithelial_cells          Fibroblasts
                   5                  729                    3
                 GMP          Hepatocytes            HSC_CD34+
                   4                  725                    4
           iPS_cells        Keratinocytes           Macrophage
                   8                    3                    4
            Monocyte                  MSC              Neurons
                   1                    6                    3
         Neutrophils              NK_cell          Osteoblasts
                   2                    7                    4
    Pre-B_cell_CD34-        Pro-Myelocyte  Smooth_muscle_cells
                   1                    1                    9
             T_cells    Tissue_stem_cells
                  26                    7

=> Epithelial_cells         
#### cluster 34 ####

Cluster34 = subset(x = immune.combined, idents = '34')

> table(Cluster34$SingleR.calls)

NK_cell T_cells
      7    1457

=> T_cells
#### cluster 35 ####

Cluster35 = subset(x = immune.combined, idents = '35')

> table(Cluster35$SingleR.calls)

           DC   Hepatocytes     HSC_CD34+    Macrophage      Monocyte
           31          1021             1           267            55
  Neutrophils       NK_cell Pro-Myelocyte       T_cells
            1             2             2            65

=> Hepatocytes    
#### cluster 36 ####

Cluster36 = subset(x = immune.combined, idents = '36')

> table(Cluster36$SingleR.calls)

          B_cell               DC              GMP      Hepatocytes
              19              284                5              106
      HSC_-G-CSF        HSC_CD34+       Macrophage         Monocyte
               2                1               17              940
         NK_cell Pre-B_cell_CD34-          T_cells
              16                4               32

=> Monocyte
#### cluster 37 ####

Cluster37 = subset(x = immune.combined, idents = '37')

> table(Cluster37$SingleR.calls)

Hepatocytes
       1221

=> Hepatocytes
#### cluster 38 ####

Cluster38 = subset(x = immune.combined, idents = '38')

> table(Cluster38$SingleR.calls)

          B_cell               DC              GMP      Hepatocytes
               1               46                1               51
      Macrophage         Monocyte          NK_cell Pre-B_cell_CD34-
             599              189               11                2
Pro-B_cell_CD34+          T_cells
               1                4

=> Macrophage        
#### cluster 39 ####

Cluster39 = subset(x = immune.combined, idents = '39')

> table(Cluster39$SingleR.calls)

          B_cell               DC              GMP      Hepatocytes
             569                1               30               12
      HSC_-G-CSF       Macrophage         Monocyte      Neutrophils
               3                2               37                1
         NK_cell Pre-B_cell_CD34- Pro-B_cell_CD34+          T_cells
              68               82                7               43

=> B_cell              
#### cluster 40 ####

Cluster40 = subset(x = immune.combined, idents = '40')

> table(Cluster40$SingleR.calls)

  Endothelial_cells         Hepatocytes             NK_cell Smooth_muscle_cells
                 36                 791                   1                   2
  Tissue_stem_cells
                  5

=> Hepatocytes            
#### cluster 41 ####

Cluster41 = subset(x = immune.combined, idents = '41')

> table(Cluster41$SingleR.calls)

             CMP      Hepatocytes       HSC_-G-CSF         Monocyte
               2                2              174               58
       Myelocyte      Neutrophils          NK_cell Pre-B_cell_CD34-
              12              472               31               22
   Pro-Myelocyte          T_cells
               3               17

=> Neutrophils         
#### cluster 42 ####

Cluster42 = subset(x = immune.combined, idents = '42')

> table(Cluster42$SingleR.calls)

        CMP          DC Hepatocytes  Macrophage    Monocyte     NK_cell
        131          12          28           2           3         106
    T_cells
        307

=> T_cells
#### cluster 43 ####

Cluster43 = subset(x = immune.combined, idents = '43')

> table(Cluster43$SingleR.calls)

Hepatocytes
        518

=> Hepatocytes
#### cluster 44 ####

Cluster44 = subset(x = immune.combined, idents = '44')

> table(Cluster44$SingleR.calls)

          CMP   Hepatocytes       NK_cell Pro-Myelocyte       T_cells
            1            14           195             1           299

=> T_cells
#### cluster 45 ####

Cluster45 = subset(x = immune.combined, idents = '45')

> table(Cluster45$SingleR.calls)

          BM   BM & Prog. Erythroblast  Hepatocytes   HSC_-G-CSF          MEP
           3            1           41           10            3            2
     NK_cell    Platelets      T_cells
          81          155          194

=> T_cells
#### cluster 46 ####

Cluster46 = subset(x = immune.combined, idents = '46')

> table(Cluster46$SingleR.calls)

Endothelial_cells       Hepatocytes
              379                42

=> Endothelial_cells      
#### cluster 47 ####

Cluster47 = subset(x = immune.combined, idents = '47')

> table(Cluster47$SingleR.calls)

          Astrocyte        Chondrocytes                  DC         Fibroblasts
                  1                   2                  17                   3
        Hepatocytes          Macrophage            Monocyte         Neutrophils
                  9                 168                   7                   1
Smooth_muscle_cells   Tissue_stem_cells
                  3                   7

=> Macrophage

```
