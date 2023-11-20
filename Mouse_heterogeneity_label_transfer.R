library(dplyr)
library(devtools)
library(Seurat)
library(Matrix)
library(AUCell)
library(GSEABase)
library(GSA)
library(DoubletFinder)
library(harmony)
library(nichenetr)
library(ggplot2)
NRAS_cell_line <- Read10X(data.dir = "/sc5rCMA309_mm10/Cellranger/Data_210415_NovaSeq2_FCB/outs/raw_feature_bc_matrix")
NRAS_cell_line <- CreateSeuratObject(counts = NRAS_cell_line, project = "NRAS_cell_line", min.cells = 10, min.features = 500)
NRAS_cell_line
NRAS_cell_line[["percent.mt"]] <- PercentageFeatureSet(NRAS_cell_line, pattern = "^mt-")
NRAS_cell_line <- SCTransform(NRAS_cell_line, verbose = TRUE, vars.to.regress = c("percent.mt"), return.only.var.genes = F)
sweep.res.list <- paramSweep_v3(NRAS_cell_line, PCs = 1:10, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pK=as.numeric(as.character(bcmvn$pK))
BCmetric=bcmvn$BCmetric
pK_choose = pK[which(BCmetric %in% max(BCmetric))]
par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
setwd("/NRAS_cell_line_B6_Nude/Results_Nude_primary")
pdf("DF_NRAS_cell_line.pdf")
plot(x = pK, y = BCmetric, pch = 16,type="b",
     col = "blue",lty=1)
abline(v=pK_choose,lwd=2,col='red',lty=2)
title("The BCmvn distributions")
text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 3,col = "red")
dev.off()
nExp_poi <- round(0.048*nrow(NRAS_cell_line@meta.data)) 
nExp_poi
NRAS_cell_line <- doubletFinder_v3(NRAS_cell_line, PCs = 1:10, pN = 0.25, pK = pK_choose, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
NRAS_cell_line$Doublets <-  NRAS_cell_line@meta.data[, grep("DF.", colnames(NRAS_cell_line@meta.data))]



NRAS_B6 <- Read10X(data.dir = "/sc5rCMA308_mm10/Cellranger/Data_210415_NovaSeq2_FCB/outs/raw_feature_bc_matrix")
NRAS_B6 <- CreateSeuratObject(counts = NRAS_B6, project = "NRAS_B6", min.cells = 10, min.features = 500)
NRAS_B6
NRAS_B6[["percent.mt"]] <- PercentageFeatureSet(NRAS_B6, pattern = "^mt-")
NRAS_B6 <- SCTransform(NRAS_B6, verbose = TRUE, vars.to.regress = c("percent.mt"), return.only.var.genes = F)
sweep.res.list <- paramSweep_v3(NRAS_B6, PCs = 1:10, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pK=as.numeric(as.character(bcmvn$pK))
BCmetric=bcmvn$BCmetric
pK_choose = pK[which(BCmetric %in% max(BCmetric))]
par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
setwd("/NRAS_cell_line_B6_Nude/Results_Nude_primary/")
pdf("DF_NRAS_B6.pdf")
plot(x = pK, y = BCmetric, pch = 16,type="b",
     col = "blue",lty=1)
abline(v=pK_choose,lwd=2,col='red',lty=2)
title("The BCmvn distributions")
text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 3,col = "red")
dev.off()
nExp_poi <- round(0.048*nrow(NRAS_B6@meta.data))  
nExp_poi
NRAS_B6 <- doubletFinder_v3(NRAS_B6, PCs = 1:10, pN = 0.25, pK = pK_choose, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
NRAS_B6$Doublets <-  NRAS_B6@meta.data[, grep("DF.", colnames(NRAS_B6@meta.data))]



NRAS_Nude  <- Read10X(data.dir = "/Volumes/Samsung_1/Endothelial_exp_invivo/TK_NRasBend3_8Tumor/raw_feature_bc_matrix")
NRAS_Nude <- CreateSeuratObject(counts = NRAS_Nude, project = "NRAS_Nude", min.cells = 10, min.features = 500)
NRAS_Nude
NRAS_Nude[["percent.mt"]] <- PercentageFeatureSet(NRAS_Nude, pattern = "^mt-")
NRAS_Nude <- SCTransform(NRAS_Nude, verbose = TRUE, vars.to.regress = c("percent.mt"), return.only.var.genes = F)
sweep.res.list <- paramSweep_v3(NRAS_Nude, PCs = 1:10, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pK=as.numeric(as.character(bcmvn$pK))
BCmetric=bcmvn$BCmetric
pK_choose = pK[which(BCmetric %in% max(BCmetric))]
par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
setwd("/NRAS_cell_line_B6_Nude/Results_Nude_primary/")
pdf("DF_NRAS_Nude.pdf")
plot(x = pK, y = BCmetric, pch = 16,type="b",
     col = "blue",lty=1)
abline(v=pK_choose,lwd=2,col='red',lty=2)
title("The BCmvn distributions")
text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 3,col = "red")
dev.off()
nExp_poi <- round(0.048*nrow(NRAS_Nude@meta.data))  
nExp_poi
NRAS_Nude <- doubletFinder_v3(NRAS_Nude, PCs = 1:10, pN = 0.25, pK = pK_choose, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
NRAS_Nude$Doublets <-  NRAS_Nude@meta.data[, grep("DF.", colnames(NRAS_Nude@meta.data))]
table(NRAS_Nude$Doublets)

All_merged <- merge(NRAS_cell_line, y=c(NRAS_B6, NRAS_Nude), project="All_merged")
All_merged
saveRDS(All_merged,"/NRAS_cell_line_B6_Nude/Results_Nude_primary/All_merged.rds")



########################################################################################################
All_merged <- readRDS("/NRAS_cell_line_B6_Nude/Results_Nude_primary/All_merged.rds")
table(All_merged$orig.ident)
DirRes <- "/NRAS_cell_line_B6_Nude/Results_Nude_primary/"

pdf(file.path(DirRes, "QC.pdf"), width = 14, height = 7)
All_merged[["percent.mt"]] <- PercentageFeatureSet(All_merged, pattern = "^mt-")
VlnPlot(All_merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0, group.by = "orig.ident")
plot1 <- FeatureScatter(All_merged, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(All_merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
pdf(file.path(DirRes, "QC_subset.pdf"), width = 14, height = 7)
All_merged <- subset(All_merged, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10)
VlnPlot( All_merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
dev.off()


mouse_cell_cycle_genes <- readRDS("/mouse_cell_cycle_genes.rds")
s.genes <- mouse_cell_cycle_genes$s.genes
g2m.genes <- mouse_cell_cycle_genes$g2m.genes
All_merged  <- CellCycleScoring(All_merged, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

############### Subset for singlets ################
All_merged_subset <- subset(All_merged, subset = Doublets == "Singlet")
All_merged_subset <- SCTransform(All_merged_subset, verbose = TRUE, vars.to.regress = c("percent.mt", 'S.Score', 'G2M.Score'))
All_merged_subset <- RunPCA(All_merged_subset , verbose = TRUE)
All_merged_subset <- RunHarmony(All_merged_subset, group.by.vars = "orig.ident", assay.use="SCT")

All_merged_subset  <- FindNeighbors(All_merged_subset , dims = 1:10, reduction = "harmony")
All_merged_subset  <- FindClusters(All_merged_subset , resolution = 0.2, reduction = "harmony")
All_merged_subset  <- RunUMAP(All_merged_subset , dims=1:10, reduction = "harmony") 
DimPlot(All_merged_subset, reduction = "umap", label = T, group.by = "seurat_clusters")
DimPlot(All_merged_subset, reduction = "umap", group.by="orig.ident")

table(All_merged_subset$seurat_clusters, All_merged_subset$orig.ident)

markers <- FindAllMarkers(All_merged_subset, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.3) # 0.3, 0.4 before
write.table(markers, "/NRAS_cell_line_B6_Nude/Results_Nude_primary/markers.txt", sep='\t', quote = FALSE, col.names = T, row.names = F)


##########################        Melanoma score         ##########################
#All_merged_subset$MLS <- NULL

genes<-c("Dct","Syt4","Ptgds","Csn3","Gpnmb","Mkln1","Mlana","4931406C07Rik","Itga4","Plp1","Lgi4","Luc7l2","Cst6","Car2","Hmga2","Kcnj10","Tyr","S100b","Atp1a1","Atp6v1f","Ube2h","Scn1b","Metrn","Cpe","Dbi","Ahcyl2","Camk2b","Tcaf1","Syngr1","Hpse","Fth1","Tecpr1","Ccnd1","Kcnma1","Akap12","Lncpint","Pla2g7","Sema3d","Kcnn4","Akr1b7","Arap2","Hip1","Slc45a2","Dhh","S100a1","Rapgef4","Cd200","Ubn2","Sox10","AC149090.1","Ccnd2","Gstp1","Etv1","C4b","Alcam","Pax3","Lgals3","Tnpo3","Plekhb1","Sort1","Ajap1","Fxyd3","Cers4","Gjc3","Tfap2b","Myo5a","Mef2c","Ndufb2","Atp1b1","Mmp16","Car6","Cavin2","Limd1","Kdm7a","Klhdc10","Erbb3","Sema3b","Serpina3n","Cdkn2b","1810058I24Rik","Tmem176a","Nrp2","Cotl1","Mest","Nav2","Mcc","Trim24","Mkrn1","Wdfy1","Mfge8","Cdh19","Shc4","Agt","Fbxl7","Lbh","Apoe","Nceh1","Mitf","Ldhb","Hipk2")
geneSets <- GeneSet(genes, setName="MLS")
geneSets
cells_rankings <- AUCell_buildRankings(All_merged_subset@assays[["SCT"]]@counts,splitByBlocks=TRUE)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)

pdf("/NRAS_cell_line_B6_Nude/Results_Nude_primary/MLS.pdf", width = 4, height = 4)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
MLS<-getAUC(cells_AUC)
MLS<-t(MLS)
All_merged_subset@meta.data<-cbind(All_merged_subset@meta.data, MLS)
FeaturePlot(All_merged_subset, features = "MLS", label = T) 
FeaturePlot(All_merged_subset, features = "nFeature_RNA", label = T) +NoAxes()
FeaturePlot(All_merged_subset, c("Pecam1", "Sox10", "Mlana", "Dct", "Lum", "Thy1", "Ptprc", "MLS"))
VlnPlot(All_merged_subset, features = "MLS") 
VlnPlot(All_merged_subset, features = "MLS", group.by =  "orig.ident") 
dev.off()
FeaturePlot(All_merged_subset, features = "Dct", label = T) +NoAxes()
VlnPlot(All_merged_subset, features = "MLS", group.by =  "orig.ident") 


hist(All_merged_subset$MLS)

################### Malignant subset
All_merged_subset_malig <- subset(All_merged_subset, subset = seurat_clusters == c("0" ,"6", "2"))
FeaturePlot(All_merged_subset_malig, features = "MLS", label = T) 
DimPlot(All_merged_subset_malig)
VlnPlot(All_merged_subset_malig, features = "MLS", group.by =  "orig.ident") 

table(All_merged_subset_malig$orig.ident)
All_merged_subset_malig <- SCTransform(All_merged_subset_malig, verbose = TRUE, vars.to.regress = c("percent.mt", 'S.Score', 'G2M.Score'))
All_merged_subset_malig <- RunPCA(All_merged_subset_malig , verbose = TRUE)
All_merged_subset_malig <- RunHarmony(All_merged_subset_malig, group.by.vars = "orig.ident", assay.use="SCT")
All_merged_subset_malig  <- FindNeighbors(All_merged_subset_malig , dims = 1:5, reduction = "harmony")
All_merged_subset_malig  <- FindClusters(All_merged_subset_malig , resolution = 0.1, reduction = "harmony")
All_merged_subset_malig  <- RunUMAP(All_merged_subset_malig , dims=1:5, reduction = "harmony") 
DimPlot(All_merged_subset_malig, reduction = "umap", label = T, group.by = "seurat_clusters")
DimPlot(All_merged_subset_malig, reduction = "umap", group.by="orig.ident")

markers <- FindAllMarkers(All_merged_subset_malig, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, group.by = "seurat_clusters")
write.table(markers, "/NRAS_cell_line_B6_Nude/Results_Nude_primary/markers_malignant.txt", sep='\t', quote = FALSE, col.names = T, row.names = F)

################# Label transfer
nras <- readRDS("/Takis_10x_NRAS_INK4A_Landscape/Malignant/NRAS13_malign_preprint.rds")
nras@meta.data$"Malignant_clusters" <- plyr::revalue(as.character(nras$seurat_clusters),
                                                                 c( "0" = "Melanocytic",
                                                                    "1" = "Neural_Crest_like",
                                                                    "2" = "RNA_processing",
                                                                    "3" = "Antigen_Presentation",
                                                                    "4" = "Stem_like",
                                                                    "5" ="Stress_hypoxia",
                                                                    "6" ="Mesenchymal_like"))

DimPlot(nras)
pancreas.features <- SelectIntegrationFeatures(object.list = c(nras, All_merged_subset_malig), nfeatures = 3000)

pancreas.anchors <- FindTransferAnchors(reference = nras, query = All_merged_subset_malig, k.anchor =20, 
                                        dims = 1:15, reduction = "pcaproject", normalization.method = "SCT", features = pancreas.features)
predictions <- TransferData(anchorset = pancreas.anchors, refdata = nras$Malignant_clusters, 
                            dims = 1:15)
All_merged_subset_malig <- AddMetaData(All_merged_subset_malig, metadata = predictions)

All_merged_subset_malig$orig.ident <-sub("_","", All_merged_subset_malig$orig.ident)
All_merged_subset_malig$orig.ident <-sub("_","", All_merged_subset_malig$orig.ident)
table(All_merged_subset_malig$orig.ident)
TEST <- All_merged_subset_malig@meta.data
cell_num <- TEST %>%
  mutate(sample_id = as.factor(paste(orig.ident,  sep="_"))) %>%
  mutate(predicted.id = as.factor(predicted.id)) %>%
  group_by(sample_id, predicted.id, .drop=FALSE) %>%
  dplyr::summarise(n=n()) %>%
  tidyr::separate(sample_id, c("orig.ident"))
cell_num
total_cells<- TEST %>%
  group_by(orig.ident) %>%
  dplyr::summarise(total = n())
total_cells
cell_percentage<- left_join(cell_num, total_cells) %>%
  mutate(percentage = n/total*100)
cell_percentage


palette.malignant <- c("#77a600","#ba5e45", "#836ab1","#b30024","grey","#20b2aa","#e02887")

pdf("/NRAS_cell_line_B6_Nude/Results_Nude_primary/bar_chart_label_transfer.pdf", width = 4, height = 6)
ggplot(cell_percentage, aes(orig.ident, y=percentage, fill=predicted.id))+geom_bar(stat="identity", position = 'stack', width = 0.9)+RotatedAxis() +  labs(fill = "NRAS_predicted_cells") + 
  theme(axis.title.x=element_blank()) + theme_classic() +
  scale_fill_manual(values=palette.malignant) +
  RotatedAxis()# +
#geom_text(aes(label = n), size = 4, hjust = 0.5, vjust = 3, position = "stack")
dev.off()
