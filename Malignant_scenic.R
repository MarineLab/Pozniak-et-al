################################################################
### Run SCENIC on scRNAseq 10X data  NextFlow & R ###
################################################################
library(Seurat)
library(SCENIC)
library(GENIE3)
library(RcisTarget)
library(AUCell)
library(loomR)
library(SCopeLoomR)
library(ComplexHeatmap)
library(ggrepel)
library(dplyr)
library(devtools)
library(Seurat)
library(Matrix)
library(AUCell)
library(GSEABase)
library(GSA)
library(dplyr)
library(MAST)
library(future)
library(ggplot2)
library(harmony)
library(DoubletFinder)
library(tidyr)
library(tibble)
library(circlize)
library(ComplexHeatmap)
library(hdf5r)

### Directories
outDirSCENIC <- "/Users/u0128760/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/MALIGNANT/SCENIC/Results"
GC_ALL_MALIGNANT <- readRDS("/Users/u0128760/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/MALIGNANT/GC_ALL_MALIGNANT.rds")
DimPlot(GC_ALL_MALIGNANT)
GC_ALL_MALIGNANT
GC_ALL_MALIGNANT_s <- GC_ALL_MALIGNANT
pyScenicLoomFile <- ("/Users/u0128760/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/MALIGNANT/SCENIC/output/GC_ALL_MALIGNANT_50x.SCENIC_SCope_output.loom")
loom <- open_loom(pyScenicLoomFile, mode="r")
loom

### problem : object as input and output don't have the same dimensions... : fixed

# Read information from loom file:
regulons_incidMat <- get_regulons(loom)
regulonsAUC <- get_regulonsAuc(loom)

regulons <- regulonsToGeneLists(regulons_incidMat)
saveRDS(regulons, "/Users/u0128760/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/MALIGNANT/SCENIC/Results/regulons.rds")
regulonsAucThresholds <- get_regulonThresholds(loom) 
embeddings <- get_embeddings(loom) 

exprMat <- get_dgem(loom)
dim(exprMat)
GC_ALL_MALIGNANT_s
cellInfo <- get_cellAnnotation(loom)
clusterings <- get_clusterings_withName(loom) 
#########################################################Heatmap 1 ##########################################
x <- as.data.frame(regulonsAUC@assays[[1]])
x1 <- as.matrix(x)
#plotTsne_AUCellApp(regulonsAUC
ht<-Heatmap(x1, 
            cluster_rows = TRUE, 
            cluster_columns = TRUE,  
            clustering_distance_columns = "euclidean",
            clustering_method_columns = "complete", 
            row_names_gp = gpar(fontsize = 2),
            show_column_names = FALSE)
pdf(paste0(outDirSCENIC,"/continous_heatmap.pdf"), height = 14, width = 7)    
draw(ht)
dev.off()

## Create binary heatmap with regulons

regulonsAucThreshold_values = as.numeric(names(regulonsAucThresholds))
BinaryRegulonsAuc = (regulonsAUC@assays[[1]] > regulonsAucThreshold_values) * 1

regulon_Thresholds <- data.frame(regulonsAucThresholds, regulonsAucThreshold_values)
row.names(regulon_Thresholds) <- NULL
write.table(regulon_Thresholds, "/Users/u0128760/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/MALIGNANT/SCENIC/Results/regulon_Thresholds.txt", sep='\t', quote = FALSE, col.names = T, row.names = F)

#########################################################Heatmap 2  ##########################################
colors = c("white", "black")
#plotTsne_AUCellApp(regulonsAUC
ht <- Heatmap(BinaryRegulonsAuc, 
              cluster_rows = TRUE, 
              cluster_columns = TRUE,  
              clustering_distance_columns = "euclidean",
              clustering_method_columns = "complete", 
              row_names_gp = gpar(fontsize = 2),
              show_column_names = FALSE,
              col = colors)
pdf(paste0(outDirSCENIC,"/binary_heatmap.pdf"), height = 14, width = 7)    
draw(ht)
dev.off()



### Integrate Seurat and SCENIC objects

## Add cluster SCENIC

clusterings$Cluster_SCENIC <- gsub("Unannotated Cluster ", "", clusterings[,1])

GC_ALL_MALIGNANT_s <- AddMetaData(
  object = GC_ALL_MALIGNANT_s,
  metadata = clusterings$Cluster_SCENIC,
  col.name = 'Cluster_SCENIC'
)


## Add regulons values

GC_ALL_MALIGNANT_s@meta.data$cellName <- rownames(GC_ALL_MALIGNANT_s@meta.data)
SCENIC_output_data <- as.data.frame(t(x))
SCENIC_output_data$cellName <- rownames(SCENIC_output_data)

colnames(SCENIC_output_data) <- gsub("_\\(\\+\\)", "_regulon", colnames(SCENIC_output_data)) # allow to distinguish regulons from TF expression when plotting features

Seurat_output_with_regulons_values <- AddMetaData(
  object = GC_ALL_MALIGNANT_s,
  metadata = SCENIC_output_data,
  col.name = colnames(SCENIC_output_data)
)
saveRDS(Seurat_output_with_regulons_values, "/Users/u0128760/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/MALIGNANT/SCENIC/Results/Seurat_output_with_regulons_values.rds")

## Add binary regulons values

SCENIC_output_data_binary <- as.data.frame(t(BinaryRegulonsAuc))
SCENIC_output_data_binary$cellName <- rownames(SCENIC_output_data_binary)

colnames(SCENIC_output_data_binary) <- gsub("_\\(\\+\\)", "_regulon", colnames(SCENIC_output_data_binary))  # allow to distinguish regulons from TF expression when plotting features

Seurat_output_with_regulons_binary_values <- AddMetaData(
  object = GC_ALL_MALIGNANT_s,
  metadata = SCENIC_output_data_binary,
  col.name = colnames(SCENIC_output_data_binary)
)
saveRDS(SCENIC_output_data_binary, "/Users/u0128760/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/MALIGNANT/SCENIC/Results/SCENIC_output_data_binary.rds")


# Add annotation to binary heatmap
Merged <- inner_join(SCENIC_output_data, GC_ALL_MALIGNANT_s@meta.data)
seurat_cluster <- Merged$seurat_clusters
scenic_cluster <- Merged$Cluster_SCENIC
orig_ident <- Merged$orig.ident


ha = HeatmapAnnotation(seurat_cluster = seurat_cluster,
                       scenic_cluster = scenic_cluster,
                       #orig_ident = orig_ident,
                       annotation_name_side = "left")

# Should add annotation to become really interesting
#           seurat_cluster = c("0" = "indianred2", "1" = "darkorange2", "2" = "yellow4",
#                              "3" = "olivedrab", "4" = "chartreuse3", "5" = "lightseagreen",
#                              "6" = "yellow", "7" = "mediumturquoise", "8" = "steelblue2",
#                              "9" = "dodgerblue2", "10" = "mediumorchid2", "11" = "mediumorchid4", 
#                              "12" = "hotpink", "13" = "deeppink2", "14" = "red"))

ht <- Heatmap(BinaryRegulonsAuc, 
              cluster_rows = TRUE, 
              cluster_columns = TRUE,  
              clustering_distance_columns = "euclidean",
              clustering_method_columns = "complete", 
              row_names_gp = gpar(fontsize = 2),
              show_column_names = FALSE,
              top_annotation = ha,
              col = colors)

pdf(paste0(outDirSCENIC,"/binary_heatmap_annotated.pdf"), height = 14, width = 7)    
draw(ht)
dev.off()


### Plot Seurat UMAP with SCENIC regulons

pdf(paste0(outDirSCENIC,"/UMAP_Seurat_with_SCENIC_Cluster.pdf"), height = 4, width = 8) 
DimPlot(GC_ALL_MALIGNANT_s, group.by = 'Cluster_SCENIC')
dev.off()

### Compare clusters of SCENIC and Seurat

pdf(paste0(outDirSCENIC,"/Comparison_cluster_Seurat_SCENIC.pdf"), height = 4, width = 6)    
ggplot(GC_ALL_MALIGNANT_s@meta.data, aes(x=seurat_clusters, fill=Cluster_SCENIC)) +
  geom_bar()
dev.off()
Seurat_output_with_regulons_values@meta.data$"Malignant_clusters" <- plyr::revalue(as.character(Seurat_output_with_regulons_values$seurat_clusters),
                                                                 c( "0" = "Melanocytic_OXPHOS",
                                                                    "1" = "Mitochondrial",
                                                                    "2" = "Melanocytic_OXPHOS",
                                                                    "3" = "Antigen_presentation",
                                                                    "4" = "Interferon_alpha_beta",
                                                                    "5" ="Trans_reg",
                                                                    "6" ="Neural_like",
                                                                    "7" ="Stress(hypoxia)",
                                                                    "8" ="Mitotic",
                                                                    "9" ="Patient_specific_A",
                                                                    "10" ="Mesenchymal",
                                                                    "11" ="Patient_specific_B"
                                                                 ))

clusters <- c(0,1,2,3,4,5,6,7,8,9,10,11)
for(i in clusters){
  list_regulons <- colnames(Seurat_output_with_regulons_values@meta.data)[196:dim(Seurat_output_with_regulons_values@meta.data)[2]]  # to check each time if right value to start with
  pvals <- c()
  ratios <- c()
  for (regulon in list_regulons) {
    # Pvals
    pval <- wilcox.test(Seurat_output_with_regulons_values@meta.data %>% dplyr::filter(seurat_clusters == i) %>% pull(regulon),
                        Seurat_output_with_regulons_values@meta.data %>% dplyr::filter(seurat_clusters != i) %>% pull(regulon))$p.value
    pvals <- c(pvals,pval)
    # Ratios
    ratio <- mean(Seurat_output_with_regulons_values@meta.data %>% dplyr::filter(seurat_clusters == i) %>% 
                    pull(regulon)) / mean(Seurat_output_with_regulons_values@meta.data %>% 
                                            dplyr::filter(seurat_clusters  != i) %>% 
                                            pull(regulon)) 
    ratios <- c(ratios,ratio)
  }
  
  pvalsAdj <- p.adjust(pvals, method = "bonferroni", n = length(pvals))
  ResultsForVolcano_0 <- as.data.frame(cbind(list_regulons,as.numeric(as.character(pvalsAdj)),as.numeric(as.character(ratios))))
  colnames(ResultsForVolcano_0) <- c('regulons','pvalsAdj','ratios')
  # Replace pvalsAdj equal to 0 by 10^-300
  ResultsForVolcano_0$pvalsAdj <- car::recode(ResultsForVolcano_0$pvalsAdj, "0 = 10^-300")
  write.table(ResultsForVolcano_0, paste0("/Users/u0128760/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/MALIGNANT/SCENIC/Results/ResultsForVolcano_", i, ".txt"), sep='\t', quote = FALSE, col.names = T, row.names = F)
}


###########################################################++++++++++++++++++++++++++++++++++++++##########
Subset_NRas <- subset(Seurat_output_with_regulons_values, subset = seurat_clusters == c("5", "6", "1"))
levels(Subset_NRas) <-  c("5", "6", "1")
levels(Subset_NRas) 

reglist <- c("seurat_clusters", "Rfx7_regulon", "Bhlha15_regulon", "Creb3l1_regulon", "Fosl1_regulon", "Gtf3a_regulon", "Msc_regulon", "Mthfd1_regulon", "Mybl1_regulon", "Nfyb_regulon", "Nkx2.9_regulon", "Sox4_regulon", "Bach1_regulon", "Phf8_regulon", "Irf1", "Irf7", "Irf9", "Nfkb2", "Stat2", "Hivep1", "Mitf", "Pax3", "Bhlhe41")
matr <- FetchData(Subset_NRas, vars = reglist, cells = NULL) # %>%
#apply(2, function(x) (x - min(x)) / (max(x) - min(x))) # %>%
#t
matr$seurat_clusters <- factor(matr$seurat_clusters, levels = c("5", "6", "1"))
matr <- matr[order(matr$seurat_clusters),]
matr$seurat_clusters
seurat_clusters <- matr$seurat_clusters
ha = HeatmapAnnotation(seurat_clusters = seurat_clusters,
                       annotation_name_side = "left")
matr$seurat_clusters <- NULL
matr <- t(matr)
mat_scaled = t(apply(matr, 1, scale))
#col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
pdf("/Users/u0128760/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/MALIGNANT/SCENIC/Results/regulons_heatmap.pdf", width = 8, height  = 7)
ht<-Heatmap(mat_scaled, 
            cluster_rows = FALSE, 
            cluster_columns = FALSE,
            show_column_names = FALSE,
            top_annotation = ha,
            #col = col_fun,
            column_split = seurat_clusters)
draw(ht)
dev.off()

VlnPlot(Seurat_output_with_regulons_values, feature=c("TCF4_regulon"), group.by = "Malignant_clusters", pt.size = 0)
FeaturePlot(Seurat_output_with_regulons_values, feature=c("MEF2C_regulon", "TFAP2B_regulon", "SOX2_regulon"), pt.size = 0)
FeaturePlot(Seurat_output_with_regulons_values, feature=c("RXRG", "TFAP2B_regulon", "SOX2_regulon"), pt.size = 0)

DimPlot(Seurat_output_with_regulons_binary_values, group.by ="_regulon")

VlnPlot(Seurat_output_with_regulons_values, feature=c("BHLHE41_regulon", "MYC_regulon", "POLR3G_regulon"), pt.size = 0)

VlnPlot(Seurat_output_with_regulons_values, feature=c("PRRX1_regulon", "MYC_regulon", "CTNNB1_regulon"), pt.size = 0)
VlnPlot(Seurat_output_with_regulons_values, feature=c("ATF4_regulon"), group.by = "orig.ident")



library(readxl)
Meta_data <- read_excel("/Users/u0128760/Documents/PROJECTS/Grand_Challenge/Meta_data/Sample_selection_Final_EL.xlsx")
Meta_data <- as.data.frame(Meta_data)
Seurat_output_with_regulons_values@meta.data$row_names <- rownames(Seurat_output_with_regulons_values@meta.data)

Meta_data <- dplyr::select(Meta_data, `Response?`, `orig.ident`)
Meta_data$Response_updated <- Meta_data$`Response?`
Meta_data$`Response?` <- NULL

dim(Seurat_output_with_regulons_values@meta.data)
Seurat_output_with_regulons_values@meta.data<-Seurat_output_with_regulons_values@meta.data %>% left_join(Meta_data, by="orig.ident")

row.names(Seurat_output_with_regulons_values@meta.data) <- Seurat_output_with_regulons_values@meta.data$x_merge
VlnPlot(Seurat_output_with_regulons_values, feature=c("TCF4_regulon"), group.by = "Malignant_clusters", split.by = "Response_updated.y")
VlnPlot(Seurat_output_with_regulons_values, feature=c("TCF4_regulon"), group.by = "Malignant_clusters", split.by = "orig.ident")
VlnPlot(Seurat_output_with_regulons_values, feature=c("TCF4"), group.by = "Malignant_clusters", split.by = "orig.ident")
