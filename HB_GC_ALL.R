## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
library(HoneyBADGER)
library(Rsamtool)
library(GenomicFeatures)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(biomaRt)
library(edgeR)
library(DESeq2)
library(ComplexHeatmap)
library(dplyr)
library(devtools)
library(Seurat)
library(Matrix)
library(AUCell)
library(dplyr)
library(circlize)
library(GSEABase)


GC_ALL_harmony_singlet <- readRDS("/GC_ALL_harmony_singlet.rds")

GC_ALL_harmony_singlet
##########################          Immune          ##########################
genes <- c("ACAP1","AKNA","ALOX5AP","ANKRD44","APOBEC3G","ARHGAP15","ARHGAP25","ARHGAP30","ARHGAP4","ARHGAP9","ARHGDIB","ATP2A3","BIN2","C16ORF54","CCDC88B","CD37","CD48","CD52","CD53","CD69","CD84","CDC42SE2","CELF2","CNTRL","CORO1A","CSK","CXCR4","CYTH4","CYTIP","DEF6","DENND1C","DOCK2","DOCK8","DUSP2","EVI2B","FERMT3","FGD3","FNBP1","GBP5","GPR65","GPSM3","HCLS1","HMHA1","IKZF1","IL10RA","IL16","IL2RG","INPP5D","ITGA4","ITGAL","ITGB2","LAIR1","LAPTM5","LCP1","LILRB3","LIMD2","LPXN","LSP1","LY9","MAP4K1","MYO1G","NCKAP1L","NR4A2","PARP8","PARVG","PIK3CD","PIM2","PLCB2","PLEKHA2","PRKCB","PSD4","PSTPIP2","PTK2B","PTPN22","PTPN6","PTPN7","PTPRC","RAC2","RASSF5","RCSD1","RGS1","RHOH","RPS6KA1","SAMSN1","SASH3","SLA","SNX20","SP140","STK17B","TAGAP","TBC1D10C","TMC6","TMC8","TMSB4X","TRAF3IP3","TSC22D3","TSTD1","UCP2","VAV1","WIPF1")
geneSets <- GeneSet(genes, setName="Immune")
geneSets
cells_rankings <- AUCell_buildRankings(GC_ALL_harmony_singlet@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Immune<-getAUC(cells_AUC)
Immune<-t(Immune)
GC_ALL_harmony_singlet@meta.data<-cbind(GC_ALL_harmony_singlet@meta.data, Immune)
#########################          Tirosh_malignant        ##########################
genes<-c("MIA","TYR","SLC45A2","CDH19","PMEL","SLC24A5","MAGEA6","GJB1","PLP1","PRAME","CAPN3","ERBB3","GPM6B","S100B","FXYD3","PAX3","S100A1","MLANA","SLC26A2","GPR143","CSPG4","SOX10","MLPH","LOXL4","PLEKHB1","RAB38","QPCT","BIRC7","MFI2","LINC00473","SEMA3B","SERPINA3","PIR","MITF","ST6GALNAC2","ROPN1B","CDH1","ABCB5","QDPR","SERPINE2","ATP1A1","ST3GAL4","CDK2","ACSL3","NT5DC3","IGSF8","MBP")
geneSets <- GeneSet(genes, setName="Lineage")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_ALL_harmony_singlet@assays[["RNA"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
pdf("/UMAP_singlets_lineage.pdf", width = 7, height = 7)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Lineage<-getAUC(cells_AUC)
Lineage<-t(Lineage)
GC_ALL_harmony_singlet@meta.data<-cbind(GC_ALL_harmony_singlet@meta.data, Lineage)
FeaturePlot(GC_ALL_harmony_singlet, features = "Tirosh_malignant")
dev.off()


################################################################# for HoneyBADGER ##############################
### annotation ###
origin <- GC_ALL_harmony_singlet@meta.data[["orig.ident"]]
immune <- GC_ALL_harmony_singlet@meta.data[["Immune"]]
malignant <- GC_ALL_harmony_singlet@meta.data[["Lineage"]]
ha = HeatmapAnnotation(origin = origin,
                       immune = immune,
                       malignant = malignant,
                       annotation_name_side = "left")

### subsetting for immune only to use as normal for HoneyBADGER ####
Immune_sub <- subset(GC_ALL_harmony_singlet, subset = Immune > 0.15)
Normal_immune <- as.matrix(GetAssayData(Immune_sub@assays[["RNA"]]))

# Load mart object for later annotation

mart.obj<-useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                  dataset = "hsapiens_gene_ensembl", 
                  host = "www.ensembl.org",
                  ensemblRedirect = FALSE)

# Load the expression data as a gene x cell matrix
exprs	<- as.matrix(GetAssayData(GC_ALL_harmony_singlet@assays[["RNA"]]))

exprs 		<- cpm(exprs, log=TRUE)

exprsNorm <- Normal_immune

exprsNorm 		<- cpm(exprsNorm, log=TRUE)

# Only retain the genes shared in both datasets
sharedGenes 	<- intersect(rownames(exprsNorm), rownames(exprs))
exprs 		<- exprs[sharedGenes,]
exprsNorm 	<-exprsNorm[sharedGenes, ]

# Keep 5000 most expressed genes (other options are possible here...)
highExpr 	<- names(head(sort(rowMeans(exprs), decreasing=TRUE), 5000))
exprs 		<- exprs[highExpr,]
exprsNorm 	<-exprsNorm[highExpr, ]
# Take average ref
exprsNorm <- rowMeans(exprsNorm)
options(future.globals.maxSize = 4000 * 1024^2)

# Run honeybadger 
hb <- new('HoneyBADGER', name='GC_BT')
hb$setGexpMats(exprs, exprsNorm, mart.obj, filter = FALSE, scale=TRUE, verbose=TRUE, id = "hgnc_symbol")
png("/HB/InitViz7500Scaled_immune.png")
hb$plotGexpProfile(gexp.norm.sub = hb$gexp.norm, setOrder = FALSE, returnPlot = FALSE)
gexpPlot <- hb$plotGexpProfile(gexp.norm.sub = hb$gexp.norm, setOrder = TRUE, returnPlot = FALSE)
dev.off()

# Run honeybadger again and create prettier plot with complexheatmap.
gexpPlot <- hb$plotGexpProfile(gexp.norm.sub = hb$gexp.norm, setOrder = FALSE, returnPlot = TRUE)
plotData <- do.call(rbind, gexpPlot)
write.table(plotData, file = "/HB/plotData_Honey_BADGER_for_score_immune.txt", quote=F, sep="\t")
saveRDS(gexpPlot, file = "/HB/HoneyBadger_gexplot.rds")

pdf("/HB/Human_SS_HoneyBadger0_2_immune.pdf")
col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
ht<-Heatmap(plotData, 
            cluster_rows = FALSE, 
            cluster_columns = TRUE,  
            clustering_distance_columns = "euclidean",
            clustering_method_columns = "complete", 
            show_column_names = FALSE,
            #heatmap_legend_param = list(legend_direction = "vertical", title_position = "leftcenter"),
            gap = unit(2, "mm"), 
            top_annotation = ha,
            col = col_fun,
            split = factor(rep(c(paste("Chr", c(1:22))), times = sapply(gexpPlot, nrow)), levels = c(paste("Chr", c(1:22)))), 
            row_title_rot = 0)
draw(ht)
dev.off()


