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
library(ComplexHeatmap)
library(circlize)
library(ggpubr)

############### read the output from CellRanger ###############
dir_list <- dir(path="/Grand_Challenge/After_Cell_Ranger/RNA", recursive = FALSE, full.names = FALSE)
setwd("/Grand_Challenge/After_Cell_Ranger/RNA")
for (i in 1:length(dir_list)) {
  setwd(file.path("/Grand_Challenge/After_Cell_Ranger/RNA", dir_list[i]))
  assign(dir_list[i], Read10X(data.dir = "outs/raw_feature_bc_matrix"))
}

# CMA135 was excluded from the analysis
############### these are names of folder outputs from CellRanger ###############
list_data <- c(scrCMA036_hg19,
               scrCMA046_hg19,
               scrCMA038_hg19,
               scrCMA044_hg19,
               scrCMA040_hg19,
               scrCMA041_hg19,
               scrCMA048_hg19,
               scrCMA049_hg19,
               scrCMA050_hg19,
               scrCMA054_hg19,
               scrCMA063_hg19,
               scrCMA055_hg19,
               scrCMA064_hg19,
               sc5rCMA061_hg19,
               sc5rCMA066_hg19,
               scrCMA068_hg19,
               scrCMA072_hg19,
               sc5rCMA070_hg19,
               sc5rCMA074_hg19,
               scrCMA076_hg19,
               scrCMA088_hg19,
               scrCMA077_hg19,
               scrCMA089_hg19,
               scrCMA087_hg19,
               scrCMA090_hg19,
               scrCMA091_hg19,
               scrCMA093_hg19,
               scrCMA112_hg19,
               scrCMA094_hg19,
               scrCMA109_hg19,
               scrCMA118_hg19,
               scrCMA119_hg19,
               scrCMA129_hg19,
               scrCMA120_hg19,
               scrCMA130_hg19,
               scrCMA121_hg19,
               scrCMA131_hg19,
               sc5rCMA136_hg19,
               sc5rCMA149_hg19,
               sc5rCMA141_hg19,
               sc5rCMA152_hg19,
               sc5rCMA155_hg19,
               sc5rCMA144_hg19,
               sc5rCMA188_hg19,
               sc5rCMA192_hg19,
               sc5rCMA196_Hg19)

############### These are my new names for seurat objects ###############
file_list1 <- c("scrCMA036",
                "scrCMA046",
                "scrCMA038",
                "scrCMA044",
                "scrCMA040",
                "scrCMA041",
                "scrCMA048",
                "scrCMA049",
                "scrCMA050",
                "scrCMA054",
                "scrCMA063",
                "scrCMA055",
                "scrCMA064",
                "sc5rCMA061",
                "sc5rCMA066",
                "scrCMA068",
                "scrCMA072",
                "sc5rCMA070",
                "sc5rCMA074",
                "scrCMA076",
                "scrCMA088",
                "scrCMA077",
                "scrCMA089",
                "scrCMA087",
                "scrCMA090",
                "scrCMA091",
                "scrCMA093",
                "scrCMA112",
                "scrCMA094",
                "scrCMA109",
                "scrCMA118",
                "scrCMA119",
                "scrCMA129",
                "scrCMA120",
                "scrCMA130",
                "scrCMA121",
                "scrCMA131",
                "sc5rCMA136",
                "sc5rCMA149",
                "sc5rCMA141",
                "sc5rCMA152",
                "sc5rCMA155",
                "sc5rCMA144",
                "sc5rCMA188",
                "sc5rCMA192",
                "sc5rCMA196")


# Create Seurat objects
for (i in 1:length(list_data)) {
  assign(file_list1[i], CreateSeuratObject(counts = list_data[[i]], project = file_list1[[i]], min.cells = 10, min.features = 1000))
}

############### This is Seurat objects list ###############
file_list2 <- c(scrCMA036,
                scrCMA046,
                scrCMA038,
                scrCMA044,
                scrCMA040,
                scrCMA041,
                scrCMA048,
                scrCMA049,
                scrCMA050,
                scrCMA054,
                scrCMA063,
                scrCMA055,
                scrCMA064,
                sc5rCMA061,
                sc5rCMA066,
                scrCMA068,
                scrCMA072,
                sc5rCMA070,
                sc5rCMA074,
                scrCMA076,
                scrCMA088,
                scrCMA077,
                scrCMA089,
                scrCMA087,
                scrCMA090,
                scrCMA091,
                scrCMA093,
                scrCMA112,
                scrCMA094,
                scrCMA109,
                scrCMA118,
                scrCMA119,
                scrCMA129,
                scrCMA120,
                scrCMA130,
                scrCMA121,
                scrCMA131,
                sc5rCMA136,
                sc5rCMA149,
                sc5rCMA141,
                sc5rCMA152,
                sc5rCMA155,
                sc5rCMA144,
                sc5rCMA188,
                sc5rCMA192,
                sc5rCMA196)


# Filter Seurat objects
for (i in 1:length(file_list2)) {
  file_list2[[i]][["percent.mt"]] <- PercentageFeatureSet(file_list2[[i]], pattern = "^MT-")
  #assign(file_list1[i], subset(file_list2[[i]], subset = nFeature_RNA > 1000 & nFeature_RNA < 7500 & percent.mt < 30))
}


# SCT normalization for all Seurat objects
for (i in 1:length(file_list2)) {
  assign(file_list1[i], SCTransform(file_list2[[i]], verbose = TRUE, vars.to.regress = c("percent.mt")))
}
list_percentages <- c(0.016,0.031,0.039,0.039,0.039,0.039,0.039,0.039,0.039,0.039,0.039,0.039,0.039,0.039,0.039,0.039,0.039,0.039,0.039,0.039,0.039,0.039,0.039,0.039,0.039,0.039,0.039,0.039,0.039,0.046,0.046,0.054,0.054,0.054,0.054,0.054,0.069,0.069,0.069,0.069,0.069,0.076,0.076,0.076,0.076,0.076,0.076)
# Doublet Finder for all Seurat objects
for (i in 1:length(file_list2)) {
  sweep.res.list <- paramSweep_v3(file_list2[[i]], PCs = 1:10, sct = TRUE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK=as.numeric(as.character(bcmvn$pK))
  BCmetric=bcmvn$BCmetric
  pK_choose = pK[which(BCmetric %in% max(BCmetric))]
  par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
  setwd("/RESULTS")
  pdf(paste(file_list1[i], ".pdf", sep=""))
  plot(x = pK, y = BCmetric, pch = 16,type="b",
       col = "blue",lty=1)
  abline(v=pK_choose,lwd=2,col='red',lty=2)
  title("The BCmvn distributions")
  text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 3,col = "red")
  dev.off()
  nExp_poi <- round(list_percentages[i]*nrow(file_list2[[i]]@meta.data))  ## Assuming list_percentages % doublet formation rate - tailor for your dataset
  nExp_poi
  assign(file_list1[i], doubletFinder_v3(file_list2[[i]], PCs = 1:10, pN = 0.25, pK = pK_choose, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE))
}

############### Update file list2 ###############
file_list2 <- c(scrCMA036,
                scrCMA046,
                scrCMA038,
                scrCMA044,
                scrCMA040,
                scrCMA041,
                scrCMA048,
                scrCMA049,
                scrCMA050,
                scrCMA054,
                scrCMA063,
                scrCMA055,
                scrCMA064,
                sc5rCMA061,
                sc5rCMA066,
                scrCMA068,
                scrCMA072,
                sc5rCMA070,
                sc5rCMA074,
                scrCMA076,
                scrCMA088,
                scrCMA077,
                scrCMA089,
                scrCMA087,
                scrCMA090,
                scrCMA091,
                scrCMA093,
                scrCMA112,
                scrCMA094,
                scrCMA109,
                scrCMA118,
                scrCMA119,
                scrCMA129,
                scrCMA120,
                scrCMA130,
                scrCMA121,
                scrCMA131,
                sc5rCMA136,
                sc5rCMA149,
                sc5rCMA141,
                sc5rCMA152,
                sc5rCMA155,
                sc5rCMA144,
                sc5rCMA188,
                sc5rCMA192,
                sc5rCMA196)
#Add Doublet finder column called Doublets in each seurat object so that when you merge the objects you have one doublets column in meta data. 
#The old DF column per sample will stay in meta data with NAs for other samples
for (i in 1:length(file_list2)) {
  file_list2[[i]]$Doublets <-  file_list2[[i]]@meta.data[, grep("DF.", colnames(file_list2[[i]]@meta.data))]
  assign(file_list1[i], file_list2[[i]])
}

GC_ALL <- merge(scrCMA036, y=c(scrCMA046,
                                        scrCMA038,
                                        scrCMA044,
                                        scrCMA040,
                                        scrCMA041,
                                        scrCMA048,
                                        scrCMA049,
                                        scrCMA050,
                                        scrCMA054,
                                        scrCMA063,
                                        scrCMA055,
                                        scrCMA064,
                                        sc5rCMA061,
                                        sc5rCMA066,
                                        scrCMA068,
                                        scrCMA072,
                                        sc5rCMA070,
                                        sc5rCMA074,
                                        scrCMA076,
                                        scrCMA088,
                                        scrCMA077,
                                        scrCMA089,
                                        scrCMA087,
                                        scrCMA090,
                                        scrCMA091,
                                        scrCMA093,
                                        scrCMA112,
                                        scrCMA094,
                                        scrCMA109,
                                        scrCMA118,
                                        scrCMA119,
                                        scrCMA129,
                                        scrCMA120,
                                        scrCMA130,
                                        scrCMA121,
                                        scrCMA131,
                                        sc5rCMA136,
                                        sc5rCMA149,
                                        sc5rCMA141,
                                        sc5rCMA152,
                                        sc5rCMA155,
                                        sc5rCMA144,
                                        sc5rCMA188,
                                        sc5rCMA192,
                                        sc5rCMA196), project="GC_ALL")

slotNames(GC_ALL)
GC_ALL
#load("~/Documents/PROJECTS/Grand_Challenge/BEFORE_TREATMENT/After_harmony.RData")
DirRes <- "/RESULTS"
GC_ALL[["percent.mt"]] <- PercentageFeatureSet(GC_ALL, pattern = "^MT-")

pdf("/RESULTS/QC.pdf", width = 7, height = 7)
VlnPlot(GC_ALL, features = c("nCount_RNA"), pt.size = 0) +  theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 5)) + NoLegend()
VlnPlot(GC_ALL, features = c("nFeature_RNA"), pt.size = 0) +  theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 5)) + NoLegend()
VlnPlot(GC_ALL, features = c("percent.mt"), pt.size = 0) +  theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 5)) + NoLegend()
FeatureScatter(GC_ALL, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(GC_ALL, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(GC_ALL, feature1 = "nFeature_RNA", feature2 = "percent.mt")
dev.off()

saveRDS(GC_ALL, "/RESULTS/GC_ALL.rds")



pdf("/RESULTS/QC_subset.pdf", width = 14, height = 7)
VlnPlot(GC_ALL_sub, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
dev.off()


give.nmedian <- function(x){
  return(c(y = median(x)*1.1, label = length(x))) 
}
pdf("/RESULTS/nCount_RNA.pdf", width = 7, height = 7)
VlnPlot(object = GC_ALL_sub, features = c("nCount_RNA"), pt.size = 0) + 
  #geom_hline(yintercept=20, linetype='dashed') +
  stat_summary(fun.data = give.nmedian, geom = "text", fun = median, size = 2) +
  theme(legend.position = "none", axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 5))
  #scale_fill_manual(values = cols)
dev.off()
#If you want it at fix position above your violin for example (ie value y=1) : 
#  give.nmedian <- function(x){
#    return(c(y = 1, label = length(x))) 
#  }
GC_ALL_harmony <- SCTransform(GC_ALL_sub, verbose = FALSE, vars.to.regress = c("percent.mt"))
GC_ALL_harmony <- RunPCA(GC_ALL_harmony , verbose = TRUE)
GC_ALL_harmony <- RunHarmony(GC_ALL_harmony, group.by.vars = "orig.ident", assay.use="SCT")

#GC_ALL_harmony <- readRDS("/RESULTS/GC_ALL_harmony.rds")
options(future.globals.maxSize = 4000 * 1024^2)
#after harmony plot harmony Embeddings on a heatmap to asses after which number the variance drops
pdf("/RESULTS/Harmony_Heatmap_all.pdf", width = 7, height = 7)
harmony_embeddings <- Embeddings(GC_ALL_harmony, 'harmony')
harmony_embeddings[1:5, 1:5]
col_fun = colorRamp2(c(-10, 0, 10), c("blue", "white", "red"))
Heatmap(harmony_embeddings, 
       cluster_rows = TRUE, 
       cluster_columns = FALSE,  
       clustering_distance_columns = "euclidean",
       clustering_method_columns = "complete", 
       show_column_names = TRUE,
       show_row_names = FALSE,
       name = "Hramony_embeedding",
       #heatmap_legend_param = list(legend_direction = "vertical", title_position = "leftcenter"),
       col = col_fun,
       row_title_rot = 0)
dev.off()

GC_ALL_harmony  <- FindNeighbors(GC_ALL_harmony , dims = 1:18, reduction = "harmony")
GC_ALL_harmony  <- FindClusters(GC_ALL_harmony , resolution = 0.75, reduction = "harmony")
GC_ALL_harmony  <- RunUMAP(GC_ALL_harmony , dims=1:18, reduction = "harmony")

# number of cells per sample
ggplot(GC_ALL_harmony_singlet@meta.data, aes(orig.ident, fill=orig.ident))+geom_bar(stat="count") + NoLegend() + theme_classic() + RotatedAxis()

pdf("/RESULTS/UMAP_all_cells_doublets.pdf", width = 7, height = 7)
DimPlot(GC_ALL_harmony, group.by = c('orig.ident'), pt.size = 0.02) +theme(legend.text = element_text(size = 5))
DimPlot(GC_ALL_harmony, group.by = c("seurat_clusters"), pt.size = 0.02, label= T) +theme(legend.text = element_text(size = 5))
DimPlot(GC_ALL_harmony, group.by = c("Doublets"), pt.size = 0.02)+theme(legend.text = element_text(size = 5))
#clculate cell cycle score and regress for it in the next step
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
GC_ALL_harmony  <- CellCycleScoring(GC_ALL_harmony, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
DimPlot(GC_ALL_harmony , reduction = "umap", group.by = 'Phase')
dev.off()
#saveRDS(GC_ALL_harmony, "/RESULTS/GC_ALL_harmony.rds")

GC_ALL_harmony <- readRDS("/RESULTS/GC_ALL_harmony.rds")

############### Subset for singlets ###############
GC_ALL_harmony
table(GC_ALL_harmony$orig.ident)
GC_ALL_harmony_singlet <- subset(GC_ALL_harmony, subset = Doublets == "Singlet")
table(GC_ALL_harmony_singlet$Doublets)
rm(GC_ALL_harmony)
############### SCT with CELL CYCLE ###############
GC_ALL_harmony_singlet
GC_ALL_harmony_singlet <- SCTransform(GC_ALL_harmony_singlet, verbose = FALSE, vars.to.regress = c("percent.mt", 'S.Score', 'G2M.Score'))
GC_ALL_harmony_singlet <- RunPCA(GC_ALL_harmony_singlet , verbose = TRUE)
GC_ALL_harmony_singlet <- RunHarmony(GC_ALL_harmony_singlet, group.by.vars = "orig.ident", assay.use="SCT")
GC_ALL_harmony_singlet

#options(future.globals.maxSize = 4000 * 1024^2)

pdf("/RESULTS/Harmony_Heatmap_all_singlets.pdf", width = 7, height = 7)
harmony_embeddings <- Embeddings(GC_ALL_harmony_singlet, 'harmony')
harmony_embeddings[1:5, 1:5]
col_fun = colorRamp2(c(-10, 0, 10), c("blue", "white", "red"))
Heatmap(harmony_embeddings, 
        cluster_rows = TRUE, 
        cluster_columns = FALSE,  
        clustering_distance_columns = "euclidean",
        clustering_method_columns = "complete", 
        show_column_names = TRUE,
        show_row_names = FALSE,
        name = "Hramony_embeedding",
        #heatmap_legend_param = list(legend_direction = "vertical", title_position = "leftcenter"),
        col = col_fun,
        row_title_rot = 0)
dev.off()

GC_ALL_harmony_singlet  <- FindNeighbors(GC_ALL_harmony_singlet , dims = 1:16, reduction = "harmony")
GC_ALL_harmony_singlet  <- FindClusters(GC_ALL_harmony_singlet , resolution = 1, reduction = "harmony")
head(Idents(GC_ALL_harmony_singlet ), 5)
GC_ALL_harmony_singlet  <- RunUMAP(GC_ALL_harmony_singlet , dims=1:16, reduction = "harmony")
pdf("/RESULTS/UMAP_singlets.pdf", width = 7, height = 7)
DimPlot(GC_ALL_harmony_singlet, group.by = c('seurat_clusters'), pt.size = 0.02, label=T) +theme(legend.text = element_text(size = 5))
DimPlot(GC_ALL_harmony_singlet, group.by = c('orig.ident'), pt.size = 0.02) +theme(legend.text = element_text(size = 5))
DimPlot(GC_ALL_harmony_singlet , reduction = "umap", group.by = 'Phase') +theme(legend.text = element_text(size = 5))
dev.off()

#saveRDS(GC_ALL_harmony_singlet, "/RESULTS/GC_ALL_harmony_singlet.rds")

GC_ALL_harmony_singlet <- readRDS("/RESULTS/GC_ALL_harmony_singlet.rds")
GC_ALL_harmony_singlet



##########################        Tirosh_malignant        ########################## should be Jerby-Arnon
genes<-c("MIA","TYR","SLC45A2","CDH19","PMEL","SLC24A5","MAGEA6","GJB1","PLP1","PRAME","CAPN3","ERBB3","GPM6B","S100B","FXYD3","PAX3","S100A1","MLANA","SLC26A2","GPR143","CSPG4","SOX10","MLPH","LOXL4","PLEKHB1","RAB38","QPCT","BIRC7","MFI2","LINC00473","SEMA3B","SERPINA3","PIR","MITF","ST6GALNAC2","ROPN1B","CDH1","ABCB5","QDPR","SERPINE2","ATP1A1","ST3GAL4","CDK2","ACSL3","NT5DC3","IGSF8","MBP")
geneSets <- GeneSet(genes, setName="Tirosh_malignant")
geneSets
cells_rankings <- AUCell_buildRankings(GC_ALL_harmony_singlet@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
pdf("/RESULTS/UMAP_singlets_tirosh_malignant.pdf", width = 7, height = 7)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Tirosh_malignant<-getAUC(cells_AUC)
Tirosh_malignant<-t(Tirosh_malignant)
GC_ALL_harmony_singlet@meta.data<-cbind(GC_ALL_harmony_singlet@meta.data, Tirosh_malignant)
FeaturePlot(GC_ALL_harmony_singlet, features = "Tirosh_malignant")
dev.off()

##########################          Immune          ##########################
genes <- c("ACAP1","AKNA","ALOX5AP","ANKRD44","APOBEC3G","ARHGAP15","ARHGAP25","ARHGAP30","ARHGAP4","ARHGAP9","ARHGDIB","ATP2A3","BIN2","C16ORF54","CCDC88B","CD37","CD48","CD52","CD53","CD69","CD84","CDC42SE2","CELF2","CNTRL","CORO1A","CSK","CXCR4","CYTH4","CYTIP","DEF6","DENND1C","DOCK2","DOCK8","DUSP2","EVI2B","FERMT3","FGD3","FNBP1","GBP5","GPR65","GPSM3","HCLS1","HMHA1","IKZF1","IL10RA","IL16","IL2RG","INPP5D","ITGA4","ITGAL","ITGB2","LAIR1","LAPTM5","LCP1","LILRB3","LIMD2","LPXN","LSP1","LY9","MAP4K1","MYO1G","NCKAP1L","NR4A2","PARP8","PARVG","PIK3CD","PIM2","PLCB2","PLEKHA2","PRKCB","PSD4","PSTPIP2","PTK2B","PTPN22","PTPN6","PTPN7","PTPRC","RAC2","RASSF5","RCSD1","RGS1","RHOH","RPS6KA1","SAMSN1","SASH3","SLA","SNX20","SP140","STK17B","TAGAP","TBC1D10C","TMC6","TMC8","TMSB4X","TRAF3IP3","TSC22D3","TSTD1","UCP2","VAV1","WIPF1")
geneSets <- GeneSet(genes, setName="Immune")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_ALL_harmony_singlet@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
pdf("/RESULTS/UMAP_singlets_immune.pdf", width = 7, height = 7)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Immune<-getAUC(cells_AUC)
Immune<-t(Immune)
GC_ALL_harmony_singlet@meta.data<-cbind(GC_ALL_harmony_singlet@meta.data, Immune)
FeaturePlot(GC_ALL_harmony_singlet , features = c ("Immune"))
dev.off()

##########################          Stroma         ##########################
genes <-c("A4GALT","ADAMTS1","ADAMTSL1","ADIRF","ANGPTL2","APP","ARHGAP29","BGN","BMP1","C1R","CCDC80","CFH","CLU","COL15A1","COL18A1","COL4A1","COL4A2","COL6A2","COX7A1","CTGF","CYB5R3","CYR61","DCHS1","DPYSL3","EFEMP1","EHD2","ELN","EPAS1","FAM171A1","FAP","FAT4","FBN1","FLRT2","FSCN1","FSTL1","GJA1","GNG11","HSPG2","HTRA1","IFITM3","IGF2","IGFBP4","IGFBP7","JAG1","KIAA1217","LAMB1","LAMB2","LAMC1","LEPROT","LHFP","LIMCH1","LIMS2","LMCD1","LOXL2","LPHN2","LRRC32","MAP1B","MEOX2","MGP","MMP2","NFIB","NID1","NNMT","NPDC1","NR2F2","NT5E","NUAK1","PEAR1","PHLDB2","PLSCR4","PPAP2A","PPAP2B","PPIC","PRKCDBP","PROCR","PRSS23","PTRF","PXDN","RAB11FIP5","RABAC1","RBPMS","RUNX1T1","S100A16","SERPINH1","SPARC","SPARCL1","STC2","TFPI","TGFB1I1","THBS1","THY1","TMEM204","TNKS1BP1","TNXB","TPBG","UNC5B","VCL","ZEB1","ZNF423","ZNF521","ICAM1","ICAM2","ICAM3","ITGA4","ITGB1","KIT","MADCAM1","MME","MMP1","MMP9","PDGFRA","PDGFRB","PECAM1","TIMP1","TIMP2","TLR1","TLR2","TLR3","TLR4","VCAM1")
geneSets <- GeneSet(genes, setName="Stroma")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_ALL_harmony_singlet@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
pdf("/RESULTS/UMAP_singlets_stroma.pdf", width = 7, height = 7)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Stroma <-getAUC(cells_AUC)
Stroma<-t(Stroma)
GC_ALL_harmony_singlet@meta.data<-cbind(GC_ALL_harmony_singlet@meta.data, Stroma )
FeaturePlot(GC_ALL_harmony_singlet , features = c ("Stroma"))
dev.off()


##########################          Endothelial_cells         ##########################
gmtFile <- paste("/Users/u0128760/Documents/PROJECTS/Grand_Challenge/Gene_signatures/Endothelial_cells.gmt")
geneSets <- getGmt(gmtFile)
geneSets
#cells_rankings <- AUCell_buildRankings(GC_ALL_harmony_singlet@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
pdf("/RESULTS/UMAP_singlets_EC.pdf", width = 7, height = 7)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Endothelial_cells <-getAUC(cells_AUC)
Endothelial_cells<-t(Endothelial_cells)
GC_ALL_harmony_singlet@meta.data<-cbind(GC_ALL_harmony_singlet@meta.data, Endothelial_cells )
FeaturePlot(GC_ALL_harmony_singlet , features = c ("Endothelial_cells"), min.cutoff = 0.052)
dev.off()


########################## CNV ##########################
#Run HB script and come back here
#Add HB results to the meta_data and calculate mean across cells
library(data.table)
cnv_honey <- fread("/HB/plotData_Honey_BADGER_for_score_immune.txt")
cnv_honey[1:5, 1:5]
dim(cnv_honey)
class(cnv_honey)
cnv_honey$V1 <- NULL
cnv_honey <- as.matrix(cnv_honey)
abs_honey <- apply(cnv_honey, 2, abs)
abs_honey[1:5, 1:5]
mean_cnv <- apply(abs_honey, 2, mean)
mean_cnv <- as.data.frame(mean_cnv)

GC_ALL_harmony_singlet@meta.data$x_merge <- rownames(GC_ALL_harmony_singlet@meta.data)
mean_cnv$x_merge <- rownames(mean_cnv)
GC_ALL_harmony_singlet@meta.data <- GC_ALL_harmony_singlet@meta.data %>% inner_join(mean_cnv, by="x_merge")
rownames(GC_ALL_harmony_singlet@meta.data) <- GC_ALL_harmony_singlet@meta.data$x_merge



########################## change clusters identity ##########################
# with keratinocytes
GC_ALL_harmony_singlet@meta.data$"Cluster2" <- plyr::revalue(as.character(GC_ALL_harmony_singlet$seurat_clusters),
                                                             c("0" = "Malignant",
                                                               "1" = "Immune",
                                                               "2" = "Malignant",
                                                               "3" = "Immune",
                                                               "4" = "Malignant_Immune_Mix",
                                                               "5" = "Immune",
                                                               "6" = "Immune",
                                                               "7" = "Stroma",
                                                               "8" = "Stroma",
                                                               "9" = "Stroma",
                                                               "10" = "Immune",
                                                               "11" = "Immune",
                                                               "12" = "Malignant",
                                                               "13" = "Immune",
                                                               "14" = "Stroma",
                                                               "15" = "Immune",
                                                               "16" = "Immune",
                                                               "17" = "Malignant",
                                                               "18" = "Immune",
                                                               "19" = "Malignant",
                                                               "20" = "Malignant",
                                                               "21" = "Keratinocytes",
                                                               "22" = "Immune",
                                                               "23" = "Immune",
                                                               "24" = "Unknown",
                                                               "25" = "Unknown",
                                                               "26" = "Immune",
                                                               "27" = "Immune",
                                                               "28" = "Malignant",
                                                               "29" = "Stroma",
                                                               "30" = "Immune"
                                                             ))
# without keratinocytes
GC_ALL_harmony_singlet@meta.data$"Cluster" <- plyr::revalue(as.character(GC_ALL_harmony_singlet$seurat_clusters),
                                                       c("0" = "Malignant",
                                                         "1" = "Immune",
                                                         "2" = "Malignant",
                                                         "3" = "Immune",
                                                         "4" = "Malignant",
                                                         "5" = "Immune",
                                                         "6" = "Immune",
                                                         "7" = "Stroma",
                                                         "8" = "Stroma",
                                                         "9" = "Stroma",
                                                         "10" = "Immune",
                                                         "11" = "Immune",
                                                         "12" = "Malignant",
                                                         "13" = "Immune",
                                                         "14" = "Stroma",
                                                         "15" = "Immune",
                                                         "16" = "Immune",
                                                         "17" = "Malignant",
                                                         "18" = "Immune",
                                                         "19" = "Malignant",
                                                         "20" = "Malignant",
                                                         "21" = "Stroma",
                                                         "22" = "Immune",
                                                         "23" = "Immune",
                                                         "24" = "Malignant",
                                                         "25" = "Malignant",
                                                         "26" = "Immune",
                                                         "27" = "Immune",
                                                         "28" = "Malignant",
                                                         "29" = "Stroma",
                                                         "30" = "Immune"
                                                       ))
DimPlot(GC_ALL_harmony_singlet, reduction = "umap", group.by = "Cluster",  cols = brewer.pal(3,"Dark2")) +NoAxes()
dev.copy2pdf(file="/RESULTS/All_cells_Annotated_clusters.pdf", width = 7.08, height = 5.8,useDingbats=FALSE,family="sans")

FeatureScatter(GC_ALL_harmony_singlet, "mean_cnv",  "Tirosh_malignant", group.by = "Cluster")
dev.copy2pdf(file="/RESULTS/scatter_mean_cnv_tirosh.pdf",useDingbats=FALSE,family="sans")

#saveRDS(GC_ALL_harmony_singlet, "/RESULTS/GC_ALL_harmony_singlet_meta_anno.rds")

# read the data again for plotting
GC_ALL_harmony_singlet <- readRDS("/RESULTS/GC_ALL_harmony_singlet_meta_anno.rds")


GC_ALL_harmony_nb_cells <- GC_ALL_harmony_singlet@meta.data %>% group_by(orig.ident) %>% summarise(nb_cells = n())
ggplot(GC_ALL_harmony_nb_cells, aes(x=factor(orig.ident, levels = orig.ident), y=nb_cells, fill = orig.ident)) + 
  geom_bar(stat = 'identity') +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10, angle=75, hjust=0.5, vjust=0.5), 
        axis.text.y = element_text(size = 16), axis.title.x = element_blank() , axis.title.y = element_text(face = "bold", size = 18) , legend.title = element_blank(), plot.title = element_text(size = 22, hjust = 0.5), legend.position='none') +
  ylab("Nb counts per cell")

ggplot(GC_ALL_harmony_singlet@meta.data, aes(x=factor(orig.ident),fill = orig.ident)) + 
  geom_bar() +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12, angle=75, hjust=0.5, vjust=0.5), axis.text.y = element_text(size = 12), axis.title.x = element_blank() , axis.title.y = element_text(face = "bold", size = 12) , legend.title = element_blank(), plot.title = element_text(size = 12, hjust = 0.5), legend.position='none') +
  #scale_fill_manual(values = cols) +
  ylab("Nb counts per cell") +
  geom_text(stat='count', aes(label=..count..), vjust=-1)


pdf("/RESULTS/nCount_ALL_per_sample.pdf", width = 14, height = 7)
give.nmedian <- function(x){
  return(c(y = median(x)*1.1, label = length(x))) 
}
VlnPlot(object = GC_ALL_harmony_singlet, features = c("nCount_RNA"), group.by = "orig.ident", pt.size = 0) + 
  #geom_hline(yintercept=20, linetype='dashed') +
  stat_summary(fun.data = give.nmedian, geom = "text", fun = median, size = 3) +
  theme(legend.position = "none", axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 5))
dev.off()
pdf("/RESULTS/nCount_ALL_per_patient.pdf", width = 14, height = 7)
give.nmedian <- function(x){
  return(c(y = median(x)*1.1, label = length(x))) 
}
VlnPlot(object = GC_ALL_harmony_singlet, features = c("nCount_RNA"), group.by = "GC number", pt.size = 0) + 
  #geom_hline(yintercept=20, linetype='dashed') +
  stat_summary(fun.data = give.nmedian, geom = "text", fun = median, size = 3) +
  theme(legend.position = "none", axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 5))
dev.off()



Idents(GC_ALL_harmony_singlet) <- "Cluster"
GC_ALL_harmony_singlet_markers <- FindAllMarkers(GC_ALL_harmony_singlet, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.3) # 0.3, 0.4 before
write.table(GC_ALL_harmony_singlet_markers, "/RESULTS/markers_stroma_immune_maligant_unknown.txt", sep='\t', quote = FALSE, col.names = T, row.names = F)


GC_ALL_harmony_singlet@meta.data$"Mutation" <- plyr::revalue(as.character(GC_ALL_harmony_singlet$`Mut type`),
                                                        c("N/A" = "NA"
                                                        ))
GC_ALL_harmony_singlet@meta.data$"Tissue" <- plyr::revalue(as.character(GC_ALL_harmony_singlet$`Tissue`),
                                                      c("Lymph node" = "LN"
                                                      ))


# read the updated meta_data
library(readxl)
Meta_data <- read_excel("/Meta_data_GEX_RESPONSE_2.0-25-08-2021.xlsx")
Meta_data <- as.data.frame(Meta_data)
GC_ALL_harmony_singlet@meta.data$row_names <- rownames(GC_ALL_harmony_singlet@meta.data)
dim(GC_ALL_harmony_singlet@meta.data)
GC_ALL_harmony_singlet@meta.data<-GC_ALL_harmony_singlet@meta.data %>% left_join(Meta_data, by="orig.ident")
row.names(GC_ALL_harmony_singlet@meta.data) <- GC_ALL_harmony_singlet@meta.data$row_names

TEST <- GC_ALL_harmony_singlet@meta.data
cell_num <- TEST %>%
  mutate(sample_id = as.factor(paste(orig.ident, `GC number`,`Tissue`,`Mutation`,`BT/OT`,`Response`,  sep="_"))) %>%
  mutate(Cluster = as.factor(Cluster)) %>%
  group_by(sample_id, Cluster, .drop=FALSE) %>%
  dplyr::summarise(n=n()) %>%
  tidyr::separate(sample_id, c("orig.ident", "GC_number", "Tissue", "Mutation", "Timepoint", "Response"))

total_cells<- TEST %>%
  group_by(orig.ident) %>%
  dplyr::summarise(total = n())
total_cells
cell_percentage<- left_join(cell_num, total_cells) %>%
  mutate(percentage = n/total*100)
cell_percentage
cell_percentage <- cell_percentage %>%
  mutate(Response = recode(Response , "0" = "NR", "1" = "R"))
cell_percentage$GC_number <- as.factor(cell_percentage$GC_number)

pdf("/All_cells_3_compartments_BT_OT.pdf", width = 17, height = 5)
ggplot(cell_percentage, aes(orig.ident, y=percentage, fill=Cluster))+geom_bar(stat="identity", position = "stack") + scale_fill_manual(values =brewer.pal(3,"Dark2")) +  facet_grid(.~GC_number+Timepoint, scale="free_x", drop = TRUE, shrink = FALSE, space = "free_x") +  theme(text=element_text(size=40)) + theme_classic()+RotatedAxis()
dev.off()
pdf("/All_cells_3_compartments_TISSUE.pdf", width = 17, height = 5)
ggplot(cell_percentage, aes(orig.ident, y=percentage, fill=Cluster))+geom_bar(stat="identity", position = "stack")+ scale_fill_manual(values =brewer.pal(3,"Dark2")) +  facet_grid(~Tissue, scale="free_x", drop = TRUE, shrink = FALSE, space = "free_x") +  theme(text=element_text(size=40)) + theme_classic()+RotatedAxis()
dev.off()
pdf("/All_cells_3_compartments_Response.pdf", width = 17, height = 5)
ggplot(cell_percentage, aes(orig.ident, y=percentage, fill=Cluster))+geom_bar(stat="identity", position = "stack") + scale_fill_manual(values =brewer.pal(3,"Dark2")) +  facet_grid(.~Response, scale="free_x", drop = TRUE, shrink = FALSE, space = "free_x") +  theme(text=element_text(size=40)) + theme_classic()+RotatedAxis()
dev.off()
