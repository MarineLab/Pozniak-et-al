library(Matrix)
library(data.table)
library(dplyr)
library(igraph)
library(RColorBrewer)
library(ggplot2)
library(patchwork)
library(dplyr)
library(AUCell)
library(GSEABase)
library(Seurat)

####Load_data#########
setwd("/Volumes/T7/VISIUM_1/Melanoma2_1/RESULTS")
data_dir <- "/Volumes/T7/VISIUM_1/Melanoma2_1/outs"
list.files(data_dir)
brain5 <-Load10X_Spatial(data.dir = data_dir)
brain5 <- PercentageFeatureSet(brain5, "^MT-", col.name = "percent_mito")
#pdf("/Volumes/T7/VISIUM_1/Melanoma2_1/RESULTS/QC_featureplot.pdf", width = 11, height = 4)
SpatialFeaturePlot(brain5, features = c("nCount_Spatial", "nFeature_Spatial", "percent_mito"))
#dev.off()
#pdf("/Volumes/T7/VISIUM_1/Melanoma2_1/RESULTS/QC_vlnplot.pdf", width = 7, height = 4)
VlnPlot(brain5, features = c("nCount_Spatial", "nFeature_Spatial", "percent_mito"), pt.size = 0.1) + NoLegend()
#dev.off()
####data_preprocessing#########
brain5 = brain5[, brain5$nFeature_Spatial > 500 & brain5$percent_mito < 3]
brain5$Spatial@data = brain5$Spatial@data[, colSums(brain5$Spatial@data != 0) > 0]

brain5<- SCTransform(brain5, assay = "Spatial", verbose = FALSE,  method = "poisson")
brain5 <- RunPCA(brain5, assay = "SCT", verbose = FALSE)
brain5 <- FindNeighbors(brain5, reduction = "pca", dims = 1:30)
brain5 <- FindClusters(brain5, verbose = FALSE)
brain5 <- RunUMAP(brain5, reduction = "pca", dims = 1:30)
SpatialDimPlot(brain5, facet.highlight = TRUE, ncol = 1,pt.size.factor = 1.0)
DefaultAssay(brain5) <- "SCT"

#brain5_markers <- FindAllMarkers(brain5, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.3)
#top10 <- brain5_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
#pdf("/Volumes/T7/VISIUM_1/Melanoma2_1/RESULTS/brain5_heatmap.pdf", width = 7, height = 10)
#DoHeatmap(brain5, features = top10$gene, group.by='seurat_clusters')
#dev.off()
#write.table(brain5_markers, "/Volumes/T7/VISIUM_1/Melanoma2_1/RESULTS/brain5_markers.txt", sep='\t', quote = FALSE, col.names = T, row.names = F)
#pdf("/Volumes/T7/VISIUM_1/Melanoma2_1/RESULTS/brain5_spatial.pdf", width = 7, height = 7)
#SpatialDimPlot(brain5, facet.highlight = TRUE, ncol = 1,pt.size.factor = 2.0)
#dev.off()

########################## prepare data for the integration
GC_all_malignant<- readRDS("/Users/u0128760/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/MALIGNANT/GC_ALL_MALIGNANT.rds")
GC_all_malignant@meta.data$"Malignant_clusters" <- plyr::revalue(as.character(GC_all_malignant$seurat_clusters),
                                                                 c( "0" = "Melanocytic",
                                                                    "1" = "Mitochondrial (low quality)",
                                                                    "2" = "Melanocytic",
                                                                    "3" = "Antigen-presenting",
                                                                    "4" = "Interferon-responsive",
                                                                    "5" ="Stress (paraspeckles)",
                                                                    "6" ="Neural Crest-like",
                                                                    "7" ="Stress (hypoxia)",
                                                                    "8" ="Mitotic",
                                                                    "9" ="Patient specific A",
                                                                    "10" ="Mesenchymal-like",
                                                                    "11" ="Patient specific B"
                                                                 ))
GC_all_malignant$Malignant_clusters <- as.factor(GC_all_malignant$Malignant_clusters)

STROMA <- readRDS("/Users/u0128760/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/STROMA/GC_stroma_FR.rds")
STROMA@meta.data$"Stromal_clust" <- plyr::revalue(as.character(STROMA$seurat_clusters),
                                                  c("0" = "Pericyte",
                                                    "1" = "CAF_POSTN",
                                                    "2" = "CAF_PLA2G2A",
                                                    "3" = "BEC", #melano
                                                    "4" = "CAF_CHI3L1",
                                                    "5" = "Pericyte",
                                                    "6" = "CAF_APCDD1",
                                                    "7" = "LEC"))
STROMA$Stromal_clust <- as.factor(STROMA$Stromal_clust)

GC_ALL_harmony_singlet <- readRDS("/Users/u0128760/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/GC_ALL_harmony_singlet_meta_anno.rds")

GC_all_immune <- readRDS("/Users/u0128760/Documents/PROJECTS/Grand_Challenge/GC_ALL_SAMPLES/RESULTS/IMMUNE/GC_ALL_immune_NEW_cluster_all_1.rds")

GC_ALL_harmony_singlet$NEW_cluster_all_1 <- GC_all_immune$NEW_cluster_all_1
#GC_ALL_harmony_singlet$NEW_cluster_all_1 <- ifelse(is.na(GC_ALL_harmony_singlet$NEW_cluster_all_1) ,GC_ALL_harmony_singlet$seurat_clusters, GC_ALL_harmony_singlet$NEW_cluster_all_1)
table(GC_ALL_harmony_singlet$NEW_cluster_all_1)
DimPlot(GC_ALL_harmony_singlet, group.by = "NEW_cluster_all_1")


GC_ALL_harmony_singlet$Detailed_cluster_all_2 <- GC_all_malignant$Malignant_clusters
GC_ALL_harmony_singlet$Detailed_cluster_all_2 <- as.character(as.factor(GC_ALL_harmony_singlet$Detailed_cluster_all_2))
DimPlot(GC_ALL_harmony_singlet, group.by = "Detailed_cluster_all_2")
GC_ALL_harmony_singlet$Detailed_cluster_all_3 <- ifelse(is.na(GC_ALL_harmony_singlet$NEW_cluster_all_1) ,GC_ALL_harmony_singlet$Detailed_cluster_all_2, GC_ALL_harmony_singlet$NEW_cluster_all_1)
DimPlot(GC_ALL_harmony_singlet, group.by = "Detailed_cluster_all_3")

GC_ALL_harmony_singlet$Detailed_cluster_all_4 <- STROMA$Stromal_clust
GC_ALL_harmony_singlet$Detailed_cluster_all_4 <- as.character(as.factor(GC_ALL_harmony_singlet$Detailed_cluster_all_4))
DimPlot(GC_ALL_harmony_singlet, group.by = "Detailed_cluster_all_4")
GC_ALL_harmony_singlet$Detailed_cluster_all_5 <- ifelse(is.na(GC_ALL_harmony_singlet$Detailed_cluster_all_4) ,GC_ALL_harmony_singlet$Detailed_cluster_all_3, GC_ALL_harmony_singlet$Detailed_cluster_all_4)
DimPlot(GC_ALL_harmony_singlet, group.by = "Detailed_cluster_all_5")


GC_ALL_harmony_singlet@meta.data$"Detailed_cluster_all_7" <- plyr::revalue(as.character(GC_ALL_harmony_singlet$Detailed_cluster_all_5),
                                                                           c("Activated_CD8_Tcells" = "CD8Tcells",
                                                                             "Antigen-presenting" = "Antigen_presenting",
                                                                             "B_cells" = "Bcells",
                                                                             "BEC" = "EC",
                                                                             "CAF_APCDD1" = "CAF",
                                                                             "CAF_CHI3L1" = "CAF",
                                                                             "CAF_PLA2G2A" = "CAF",
                                                                             "CAF_POSTN" = "CAF",
                                                                             "CD4_Tcells" = "CD4Tcells",
                                                                             "Cycling_CD8_Tcells" = "NA",
                                                                             "Cytotoxic_CD8_Tcells" = "CD8Tcells",
                                                                             "DC1" = "DCs",
                                                                             "DC2" = "DCs",
                                                                             "Doublets" = "NA",
                                                                             "Dysfuntional_CD8_Tcells" = "CD8Tcells",
                                                                             "Interferon-responsive" = "Interferon_responsive",
                                                                             "LEC" = "EC",
                                                                             "M2 Macrophages" = "Macrophages",
                                                                             "Macrophages_CXCL10" = "Macrophages",
                                                                             "Macrophages_FOS" = "Macrophages",
                                                                             "Macrophages_necrosis" = "Macrophages",
                                                                             "Macrophages_PLA2G2D" = "Macrophages",
                                                                             "Mast_cells" = "NA",
                                                                             "Melanocytic" = "Melanocytic",
                                                                             "Melanoma_immune_like" = "NA",
                                                                             "Melanophages" = "NA",
                                                                             "Memory_CD8_Tcells" = "CD8Tcells",
                                                                             "Mesenchymal-like" = "Mesenchymal_like",
                                                                             "Mitochondrial (low quality)" = "NA",
                                                                             "Mitotic" = "Mitotic",
                                                                             "Monocytes_CD14" = "Monocytes",
                                                                             "Monocytes_CD16" = "Monocytes",
                                                                             "Monocytes_TNFSF13" = "Monocytes",
                                                                             "Neural Crest-like" = "Neural_Crest_like",
                                                                             "NK" = "NK",
                                                                             "Patient specific A" = "NA",
                                                                             "Patient specific B" = "NA",
                                                                             "pDC" = "pDC",
                                                                             "Pericyte" = "EC",
                                                                             "Plasma_cells" = "PlasmaCells",
                                                                             "Stress (hypoxia)" = "Stress_hypoxia",
                                                                             "Stress (paraspeckles)" = "Stress_paraspeckles",
                                                                             "Tregs" = "Tregs"))

allen_reference2 <- GC_ALL_harmony_singlet
table(allen_reference2$Detailed_cluster_all_7)
allen_reference2 <- subset(allen_reference2, subset = Detailed_cluster_all_7 !="NA")
allen_reference2 <- subset(allen_reference2, subset = `BT/OT` !="OT")
allen_reference2 <- subset(allen_reference2, subset = `Tissue` =="Lymph node")

table(allen_reference2$Tissue)

allen_reference2 <- SCTransform(allen_reference2, verbose = FALSE, vars.to.regress = c("percent.mt", 'S.Score', 'G2M.Score'), variable.features.n = 3000)
#allen_reference2$Detailed_cluster_all_7 <- as.factor(as.character(allen_reference2$Detailed_cluster_all_7))
DefaultAssay(brain5) <- "SCT"

anchors <- FindTransferAnchors(reference = allen_reference2, query = brain5, normalization.method = c("SCT"), reduction = "cca",  k.anchor = 10)
anchors
anchors_extract <- as.data.frame(anchors@anchor.features)
predictions.assay <- TransferData(anchorset = anchors, refdata = allen_reference2$Detailed_cluster_all_7, prediction.assay = TRUE, weight.reduction = "cca", dims = 1:30)
brain5[["predictions"]] <- predictions.assay
DefaultAssay(brain5) <- "predictions"
pdf("/Volumes/T7/VISIUM_1/Melanoma2_1/RESULTS/hig_res_brain5_predicted_id_melanoma.pdf", width = 15, height = 10)
SpatialFeaturePlot(brain5, features = c("Antigen-presenting","Interferon-responsive","Melanocytic","Mesenchymal-like","Neural-Crest-like","Stress-hypoxia", "Stress-paraspeckles"), 
                   pt.size.factor = 2.0, ncol = 4, crop = TRUE)
dev.off()

pdf("/Volumes/T7/VISIUM_1/Melanoma2_1/RESULTS/hig_res_brain5_predicted_id_tme.pdf", width = 15, height = 10)
SpatialFeaturePlot(brain5, features = c("CD8Tcells","CD4Tcells", "Tregs", "NK", "Bcells", "Macrophages", "Monocytes","DCs","PlasmaCells", "EC", "CAF"), 
                   pt.size.factor = 2.0, ncol = 4, crop = TRUE)
dev.off()

hist(brain5$predictions@data[1,])
##########################################################################################
################## from here CELLTREK
allen_reference2@meta.data <- allen_reference2@meta.data %>% dplyr::select(Detailed_cluster_all_7, orig.ident)

brain5 <- RenameCells(brain5, new.names=make.names(Cells(brain5)))
allen_reference2 <- RenameCells(allen_reference2, new.names=make.names(Cells(allen_reference2)))

## Visualize the ST data
SpatialDimPlot(brain5)

#allen_reference2$Detailed_cluster_all_7 <- as.factor(allen_reference2$Detailed_cluster_all_7)
brain_traint <- CellTrek::traint(st_data=brain5, sc_data=allen_reference2, sc_assay='RNA', cell_names='Detailed_cluster_all_7')
DimPlot(brain_traint, group.by = "type") 
DimPlot(brain_traint, group.by = "Detailed_cluster_all_7") 
table(brain_traint$Detailed_cluster_all_7)
list(brain_traint$Detailed_cluster_all_7)

brain5_melanoma_sub_celltrek <- CellTrek::celltrek(st_sc_int=brain_traint, int_assay='traint', sc_data=allen_reference2, sc_assay = 'RNA', 
                                                   reduction='pca', intp=T, intp_pnt=5000, intp_lin=F, nPCs=30, ntree=1000, 
                                                   dist_thresh=0.55, top_spot=5, spot_n=5, repel_r=30, repel_iter=20, keep_model=F)$celltrek

#CellTrek::celltrek_vis(brain5_melanoma_sub_celltrek@meta.data %>% dplyr::select(coord_x, coord_y, Detailed_cluster_all_7:id_new),
 #                                                             brain5_melanoma_sub_celltrek@images$slice1@image, brain5_melanoma_sub_celltrek@images$slice@scale.factors$lowres)

glut_cell <- c( "Antigen_presenting",
                "Bcells",
                "CAF",
                "CD4Tcells",
                "CD8Tcells",
                "DCs",
                "EC",
                "Interferon_responsive",
                "Macrophages",
                "Melanocytic",
                "Mesenchymal_like",
                "Mitotic",
                "Monocytes",
                "Neural_Crest_like",
                "NK",
                "pDC",
                "PlasmaCells",
                "Stress_hypoxia",
                "Stress_paraspeckles",
                "Tregs")
names(glut_cell) <- make.names(glut_cell)
brain_celltrek_glut <- brain5_melanoma_sub_celltrek
brain_celltrek_glut$Detailed_cluster_all_7 <- factor(brain_celltrek_glut$Detailed_cluster_all_7)
brain_sgraph_KL <- CellTrek::scoloc(brain_celltrek_glut, col_cell='Detailed_cluster_all_7', use_method='DT')

brain_sgraph_KL_mst_cons <- brain_sgraph_KL$mst_cons
rownames(brain_sgraph_KL_mst_cons) <- colnames(brain_sgraph_KL_mst_cons) <- glut_cell[colnames(brain_sgraph_KL_mst_cons)]
## We then extract the metadata (including cell types and their frequencies)
brain_cell_class <- brain5_melanoma_sub_celltrek@meta.data %>% dplyr::select(id=Detailed_cluster_all_7) %>% unique
brain_celltrek_count <- data.frame(freq = table(brain5_melanoma_sub_celltrek$Detailed_cluster_all_7))
brain_cell_class_new <- merge(brain_cell_class, brain_celltrek_count, by.x ="id", by.y = "freq.Var1")

CellTrek::scoloc_vis(brain_sgraph_KL_mst_cons, meta_data=brain_cell_class_new)





save.image("/Volumes/T7/VISIUM_1/Melanoma2_1/RESULTS/Malignant/High_Res_after_cellTREK.RData")




#load("/Volumes/T7/VISIUM_1/Melanoma2_1/RESULTS/Malignant/High_Res_after_cellTREK.RData")
path_to_results <- ("/Volumes/T7/VISIUM_1/Melanoma2_1/RESULTS/HIGH_RES/")
#CellTrek::scoloc_vis(brain_sgraph_KL_mst_cons, meta_data=brain_cell_class_new)

colourCount = length(unique(brain5_melanoma_sub_celltrek$Detailed_cluster_all_7))
mycolors <- colorRampPalette(brewer.pal(9, "Set1"))(colourCount)
                                                                                                                                                      
                                                                                                                                                    
pdf(paste0(path_to_results, "Spatial_plot.pdf"))
SpatialDimPlot(brain5_melanoma_sub_celltrek, group.by = "Detailed_cluster_all_7", image.alpha = 0.7, pt.size.factor = 1.8)+scale_fill_manual(values = mycolors)
dev.off()


test1<- CellTrek::run_kdist(brain_celltrek_glut, grp_col='Detailed_cluster_all_7',ref=c( "Antigen_presenting",
                                                                                         "Bcells",
                                                                                         "CAF",
                                                                                         "CD4Tcells",
                                                                                         "CD8Tcells",
                                                                                         "DCs",
                                                                                         "EC",
                                                                                         "Interferon_responsive",
                                                                                         "Macrophages",
                                                                                         "Melanocytic",
                                                                                         "Mesenchymal_like",
                                                                                         "Mitotic",
                                                                                         "Monocytes",
                                                                                         "Neural_Crest_like",
                                                                                         "NK",
                                                                                         "pDC",
                                                                                         "PlasmaCells",
                                                                                         "Stress_hypoxia",
                                                                                         "Stress_paraspeckles",
                                                                                         "Tregs"), ref_type='each',que=unique(brain_celltrek_glut$Detailed_cluster_all_7), keep_nn=T)
str(test1)



test1 <- test1 %>% subset(Detailed_cluster_all_7 %in% c("Antigen_presenting",
                                                        "Interferon_responsive",
                                                        "Melanocytic",
                                                        "Mesenchymal_like",
                                                        "Mitotic",
                                                        "Neural_Crest_like",
                                                        "Stress_hypoxia",
                                                        "Stress_paraspeckles"))

distances_list2 <- c(  "Antigen_presenting_kdist",
                       "Bcells_kdist",
                       "CAF_kdist",
                       "CD4Tcells_kdist",
                       "CD8Tcells_kdist",
                       "DCs_kdist",
                       "EC_kdist",
                       "Interferon_responsive_kdist",
                       "Macrophages_kdist",
                       "Melanocytic_kdist",
                       "Mesenchymal_like_kdist",
                       "Mitotic_kdist",
                       "Monocytes_kdist",
                       "Neural_Crest_like_kdist",
                       "NK_kdist",
                       "pDC_kdist",
                       "PlasmaCells_kdist",
                       "Stress_hypoxia_kdist",
                       "Stress_paraspeckles_kdist",
                       "Tregs_kdist")
for (i in distances_list2){
  test1[[i]] <- test1[[i]] / max(test1[[i]])
}

for (i in distances_list2) {
  
  xplot<-  VlnPlot(test1, i, group.by = 'Detailed_cluster_all_7', pt.size = 0, sort = TRUE)+  geom_boxplot(width=0.1, fill="white", outlier.size = 0, notch = T)+NoLegend()
  pdf(paste0(path_to_results, i, ".pdf"))
  print(xplot)
  dev.off()
}

################################## HEATMAP of K-distance
distances_list2 <- c("Detailed_cluster_all_7", "Antigen_presenting_kdist",
                     "Bcells_kdist",
                     "CAF_kdist",
                     "CD4Tcells_kdist",
                     "CD8Tcells_kdist",
                     "DCs_kdist",
                     "EC_kdist",
                     "Interferon_responsive_kdist",
                     "Macrophages_kdist",
                     "Melanocytic_kdist",
                     "Mesenchymal_like_kdist",
                     "Mitotic_kdist",
                     "Monocytes_kdist",
                     "Neural_Crest_like_kdist",
                     "NK_kdist",
                     "pDC_kdist",
                     "PlasmaCells_kdist",
                     "Stress_hypoxia_kdist",
                     "Stress_paraspeckles_kdist",
                     "Tregs_kdist")
matr <- FetchData(test1, vars = distances_list2, cells = NULL) # %>%
#apply(2, function(x) (x - min(x)) / (max(x) - min(x))) # %>%
#t
library(scales)
matr <- matr[ , !duplicated(colnames(matr))]
matr <- matr %>% group_by(Detailed_cluster_all_7) %>% summarise_all(mean)
#


Detailed_cluster_all_7 <- levels(test1$Detailed_cluster_all_7)
Detailed_cluster_all_7
#color_1 <- setNames(hue_pal()(11), Detailed_cluster_all_7)
color_list <- list(color_1 = Detailed_cluster_all_7) 
color_list
#names(color_list) <- list(Malignant_clusters_old)
matr$Detailed_cluster_all_7 <- factor(matr$Detailed_cluster_all_7, levels=c("Antigen_presenting",
                                                                            "Bcells",
                                                                            "CAF",
                                                                            "CD4Tcells",
                                                                            "CD8Tcells",
                                                                            "DCs",
                                                                            "EC",
                                                                            "Interferon_responsive",
                                                                            "Macrophages",
                                                                            "Melanocytic",
                                                                            "Mesenchymal_like",
                                                                            "Mitotic",
                                                                            "Monocytes",
                                                                            "Neural_Crest_like",
                                                                            "NK",
                                                                            "pDC",
                                                                            "PlasmaCells",
                                                                            "Stress_hypoxia",
                                                                            "Stress_paraspeckles",
                                                                            "Tregs"))
Detailed_cluster_all_7 <- matr$Detailed_cluster_all_7
ha = ComplexHeatmap::HeatmapAnnotation(Detailed_cluster_all_7 = Detailed_cluster_all_7,
                                       annotation_name_side = "left",
                                       annotation_name_rot = 180,
                                       show_annotation_name = FALSE
                                       #col = color_list
)



head(matr)

matr$Detailed_cluster_all_7<- NULL
head(matr)

#mat_num <- matrix(as.numeric(as.character((matr),    # Convert to numeric matrix
#                                         ncol = ncol(matr))))
matr <- t(matr)
#mat_scaled = t(apply(matr, 1, scale))
ht<-ComplexHeatmap::Heatmap(matr, 
                            cluster_rows = FALSE, 
                            cluster_columns = FALSE,
                            show_column_names = FALSE,
                            top_annotation = ha,
                            #col = color_list,
                            column_split = Detailed_cluster_all_7,
                            column_title_rot = 90)
pdf(paste0(path_to_results, "HEATMAP_distance.pdf"), width=10, height = 10)
print(ht)
dev.off()