library(dplyr)
library(devtools)
library(Seurat)
library(Matrix)
library(AUCell)
library(GSEABase)
library(GSA)
library(hypeR)
library(grkmisc)
library(ggpubr)
library(purrr)
library(ggrepel)
library(plyr)

GC_ALL_harmony_singlet <- readRDS("GC_ALL_harmony_singlet_meta_anno.rds")
DimPlot(GC_ALL_harmony_singlet, label = T)

#Identification of Melanoma Score genes
SuperMEL <- FindMarkers(GC_ALL_harmony_singlet, ident.1 = c(2,28,0,12,17,20,19) , ident.2 = c(7,8),  only.pos = T,  min.pct = 0.3, logfc.threshold = 0.3)
write.table(SuperMEL, "SuperMEL.txt", sep='\t', quote = FALSE, row.names=T, col.names = T)

######################### Melanoma_specific_genes_vs_CAFs   - "Melanoma Score" in the paper ##########################
VlnPlot(GC_ALL_harmony_singlet, feature = c("EDNRB", "MYO10","PLP1","ERBB3", "SYNGR1")) #specific only to melanoma cells
#specific only to melanoma cells
genes<-c("EDNRB", "MYO10","PLP1", "ERBB3","SYNGR1")
geneSets <- GeneSet(genes, setName="Melanoma_specific_genes_vs_CAFs")
geneSets
cells_rankings <- AUCell_buildRankings(GC_ALL_harmony_singlet@assays[["RNA"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
pdf("MALIGNANT/UMAP_Melanoma_specific_genes_vs_CAFs.pdf", width = 7, height = 7)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Melanoma_specific_genes_vs_CAFs<-getAUC(cells_AUC)
Melanoma_specific_genes_vs_CAFs<-t(Melanoma_specific_genes_vs_CAFs)
GC_ALL_harmony_singlet@meta.data<-cbind(GC_ALL_harmony_singlet@meta.data, Melanoma_specific_genes_vs_CAFs)
FeaturePlot(GC_ALL_harmony_singlet, features = "Melanoma_specific_genes_vs_CAFs")
dev.off()


######################### Subset the cells for only malignant ones #########################
GC_all_malignant <- subset(GC_ALL_harmony_singlet, subset = Tirosh_malignant > 0.11 | mean_cnv > 0.15)
GC_all_malignant <- subset(GC_all_malignant, subset = Melanoma_specific_genes_vs_CAFs > 0.2 | mean_cnv > 0.15)
GC_all_malignant <- subset(GC_all_malignant, subset = PTPRC < 0.0001)
GC_all_malignant
DimPlot(GC_all_malignant)
ggplot(GC_all_malignant@meta.data, aes(`GC number`))+geom_bar(stat="count")+facet_wrap(~GC_all_malignant$Response + GC_all_malignant@meta.data$`BT/OT`)
table(GC_all_malignant$orig.ident)

#remove samples with less than 10 cells
library(purrr)
`%not_in%` <- negate(`%in%`) 
GC_all_malignant <- GC_all_malignant %>% subset(orig.ident  %not_in% c("sc5rCMA066", "scrCMA109", "scrCMA112", "scrCMA118"))
GC_all_malignant

table(GC_all_malignant$orig.ident)
table(GC_ALL_harmony_singlet$orig.ident)

GC_all_malignant <- SCTransform(GC_all_malignant, verbose = FALSE, vars.to.regress = c("percent.mt", 'S.Score', 'G2M.Score'))
GC_all_malignant <- RunPCA(GC_all_malignant , verbose = TRUE)
GC_all_malignant <- RunHarmony(GC_all_malignant, group.by.vars = "orig.ident", assay.use="SCT")
#saveRDS(GC_all_malignant, "MALIGNANT/GC_ALL_MALIGNANT_interm.rds")
#GC_all_malignant <- readRDS("MALIGNANT/GC_ALL_MALIGNANT_interm.rds")
options(future.globals.maxSize = 4000 * 1024^2)
pdf("MALIGNANT/Harmony_Heatmap_Malignant.pdf", width = 7, height = 7)
harmony_embeddings <- Embeddings(GC_all_malignant, 'harmony')
harmony_embeddings[1:5, 1:5]
col_fun = colorRamp2(c(-10, 0, 10), c("blue", "white", "red"))
ComplexHeatmap::Heatmap(harmony_embeddings, 
    cluster_rows = TRUE, 
    cluster_columns = FALSE,  
    clustering_distance_columns = "euclidean",
    clustering_method_columns = "complete", 
    show_column_names = TRUE,
    show_row_names = FALSE,
    name = "Hramony_embeedding",
    #heatmap_legend_param = list(legend_direction = "vertical", title_position = "leftcenter"),
    #col = col_fun,
    row_title_rot = 0)
dev.off()

GC_all_malignant  <- FindNeighbors(GC_all_malignant , dims = 1:17, reduction = "harmony")
GC_all_malignant  <- FindClusters(GC_all_malignant , resolution = 0.35, reduction = "harmony")
GC_all_malignant  <- RunUMAP(GC_all_malignant , dims=1:17, reduction = "harmony")
DimPlot(GC_all_malignant, group.by = c('seurat_clusters'), label = T)

pdf("MALIGNANT/UMAP_clusters.pdf", width = 7.08, height = 5.8)
DimPlot(GC_all_malignant, group.by = c('seurat_clusters'), label = T) +NoAxes()
dev.off()


pdf("MALIGNANT/UMAP_BT_OT.pdf", width = 4, height = 3)
DimPlot(GC_all_malignant, group.by = "BT/OT", cols = c("#00AFBB", "#FC4E07"))  & NoAxes()
dev.off()

#saveRDS(GC_all_malignant, "MALIGNANT/GC_ALL_MALIGNANT.rds")

####################################### Cluster Identity #########################
GC_all_malignant_markers <- FindAllMarkers(GC_all_malignant, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.3)
write.table(GC_all_malignant_markers, "MALIGNANT/GC_all_malignant_markers.txt", sep='\t', quote = FALSE, col.names = T, row.names = F)
top20 <- GC_all_malignant_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)

pdf("MALIGNANT/Heatmap_GC_all_malignant.pdf", width = 15, height = 25)
DoHeatmap(GC_all_malignant, features = top20$gene) + NoLegend()
dev.off()

####################################### STATS - malignant clusters vs Response #######################################
GC_all_malignant<- readRDS("MALIGNANT/GC_ALL_MALIGNANT.rds")

#Add identity of the clusters
GC_all_malignant@meta.data$"Malignant_clusters" <- plyr::revalue(as.character(GC_all_malignant$seurat_clusters),
                                                                c( "0" = "Melanocytic",
                                                                   "1" = "Mitochondrial(low_quality)",
                                                                  "2" = "Melanocytic",
                                                                  "3" = "Antigen_Presentation",
                                                                   "4" = "Interferon_Alpha_Beta_Response",
                                                                   "5" ="Stress (p53 Response)",
                                                                   "6" ="Neural_Crest_like",
                                                                    "7" ="Stress (Hypoxia Response)",
                                                                   "8" ="Mitotic",
                                                                   "9" ="Patient_specific_A",
                                                                   "10" ="Mesenchymal_like",
                                                                   "11" ="Patient_specific_B"
                                                                ))

#Idents(GC_all_malignant) <- "Malignant_clusters"
#GC_all_malignant_markers <- FindAllMarkers(GC_all_malignant, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.3) # 0.3, 0.4 before
#write.table(GC_all_malignant_markers, "MALIGNANT/GC_all_malignant_markers_annotated_cluster.txt", sep='\t', quote = FALSE, col.names = T, row.names = F)


# add information about TILs infitration pattern
GC_all_malignant@meta.data$"TILs" <- plyr::revalue(as.character(GC_all_malignant$orig.ident),
                                      c("scrCMA036"="NA",
                                        "scrCMA046" = "NA",
                                        "scrCMA038" = "NonBrisk",
                                        "scrCMA044" = "Brisk",
                                        "scrCMA040" = "Absent",
                                        "scrCMA048" = "NonBrisk",
                                        "scrCMA049" = "NonBrisk",
                                        "scrCMA050" = "NonBrisk",
                                        "scrCMA054" = "Absent",
                                        "scrCMA063" = "Absent",
                                        "scrCMA055" = "Absent",
                                        "scrCMA064" = "NonBrisk",
                                        "sc5rCMA061" = "NA",
                                        "sc5rCMA066" = "Absent",
                                        "scrCMA068" = "NonBrisk",
                                        "scrCMA072" = "NonBrisk",
                                        "sc5rCMA070" = "NonBrisk",
                                        "sc5rCMA074" = "NA",
                                        "scrCMA076" = "NonBrisk",
                                        "scrCMA088" = "Brisk",
                                        "scrCMA077" = "NonBrisk",
                                        "scrCMA087" = "NA",
                                        "scrCMA090" = "Absent",
                                        "scrCMA091" = "NonBrisk",
                                        "scrCMA093" = "Necrosis",
                                        "scrCMA112" = "Necrosis",
                                        "scrCMA094" = "NA",
                                        "scrCMA109" = "NoTumor",
                                        "scrCMA119" = "Brisk",
                                        "scrCMA129" = "Necrosis",
                                        "scrCMA120" = "Absent",
                                        "scrCMA121" = "NonBrisk",
                                        "scrCMA131" = "Brisk",
                                        "sc5rCMA136" = "NonBrisk",
                                        "sc5rCMA149" = "NonBrisk",
                                        "sc5rCMA141" = "NonBrisk",
                                        "sc5rCMA152" = "NA",
                                        "sc5rCMA144" = "Brisk",
                                        "sc5rCMA155" = "NA",
                                        "sc5rCMA188" = "NA",
                                        "scrCMA089" = "NA",
                                        "scrCMA118" = "NA",
                                        "scrCMA130" = "NA",
                                        "scrCMA041" = "NA"
                                      ))
table(GC_all_malignant@meta.data$"TILs" )
###################### percentages  out of malignant cells ##################
###########################################################################

GC_all_malignant@meta.data$"TILs_bin" <- plyr::revalue(as.character(GC_all_malignant$TILs),
                                                                c("Absent" = "Absent",
                                                                  "NonBrisk" = "Present",
                                                                  "Brisk" = "Present"
                                                                ))
GC_all_malignant@meta.data$"Mutation" <- plyr::revalue(as.character(GC_all_malignant$`Mut type`),
                                                                c("N/A" = "NA"
                                                                ))
GC_all_malignant@meta.data$"Tissue" <- plyr::revalue(as.character(GC_all_malignant$`Tissue`),
                                                              c("Lymph node" = "LN"
                                                              ))
GC_all_malignant@meta.data$"Treatment" <- plyr::revalue(as.character(GC_all_malignant$`Treatment`),
                                                     c("Nivolumab + Ipilimumab" = "NNivolumabIpilimumab"
                                                     ))
############### meta data updated response
library(readxl)
Meta_data <- read_excel("/Meta_data_GEX_RESPONSE_2.0-25-08-2021.xlsx")
Meta_data <- as.data.frame(Meta_data)
GC_all_malignant@meta.data$row_names <- rownames(GC_all_malignant@meta.data)
dim(GC_all_malignant@meta.data)
GC_all_malignant@meta.data<-GC_all_malignant@meta.data %>% left_join(Meta_data, by="orig.ident")
row.names(GC_all_malignant@meta.data) <- GC_all_malignant@meta.data$row_names
table(GC_all_malignant$orig.ident)

TEST <- GC_all_malignant@meta.data
cell_num <- TEST %>%
  mutate(sample_id = as.factor(paste(orig.ident,`Treatment`, `Response`,`BT/OT`, `GC number`, `TILs_bin`, `TILs`, `Tissue`, `Mutation`,  sep="_"))) %>%
  mutate(Malignant_clusters = as.factor(Malignant_clusters)) %>%
  group_by(sample_id, Malignant_clusters, .drop=FALSE) %>%
  dplyr::summarise(n=n()) %>%
  tidyr::separate(sample_id, c("orig.ident", "Treatment", "Response", "Timepoint", "GC_number", "TILs_bin", "TILs", "Tissue", "Mutation_type"))

total_cells<- TEST %>%
  group_by(orig.ident) %>%
  dplyr::summarise(total = n())
total_cells
cell_percentage<- left_join(cell_num, total_cells) %>%
  mutate(percentage = n/total*100)
cell_percentage

cell_percentage <- subset(cell_percentage, Response != "NA")

cell_percentage <- cell_percentage %>%
  mutate(Response = recode(Response , "0" = "NR", "1" = "R"))

`%not_in%` <- negate(`%in%`) 

#cell_percentage_no <-  cell_percentage %>% subset(Malignant_clusters %not_in% c("Antigen_presenting", "Mesenchymal_like"))
cell_percentage_no <-  cell_percentage %>% subset(Malignant_clusters %not_in% c("Patient_specific_A", "Patient_specific_B"))

ggboxplot(cell_percentage_no, x = "Response", y = "percentage",
          fill = "Response", palette =c( "#FC4E07", "#00AFBB"),
          shape = "Response",
          add = "jitter")+
  stat_compare_means(aes(group = Response, label = paste0("p = ", ..p.format..), method =  "wilcox.test"),label.x = c(1))  +
  facet_wrap(~Malignant_clusters, nrow = 2, scales = "free_y") + theme_classic() + theme(text=element_text(size=10))
dev.copy2pdf(file="MALIGNANT/Cell_percentages_RvsNR.pdf", width = 15, height = 6,useDingbats=FALSE,family="sans")

cell_percentage_A <- cell_percentage %>% subset(Malignant_clusters %in% c("Antigen_Presentation"))
ggboxplot(cell_percentage_A, x = "Response", y = "percentage",
          fill = "Response", palette =c( "#FC4E07", "#00AFBB"),
          shape = "Response",
          add = "jitter")+
  stat_compare_means(aes(group = Response, label = paste0("p = ", ..p.format..), method =  "wilcox.test"),label.x = c(1), label.y = c(74))  +
  facet_wrap(~Malignant_clusters, scales = "free_y") + theme_classic() + theme(text=element_text(size=15)) + NoLegend()
dev.copy2pdf(file="MALIGNANT/Antigen_Cell_percentages_RvsNR.pdf", width = 4, height = 5,useDingbats=FALSE,family="sans")

cell_percentage_M <- cell_percentage %>% subset(Malignant_clusters %in% c("Mesenchymal_like"))
ggboxplot(cell_percentage_M, x = "Response", y = "percentage",
          fill = "Response", palette =c( "#FC4E07", "#00AFBB"),
          shape = "Response",
          add = "jitter")+
  stat_compare_means(aes(group = Response, label = paste0("p = ", ..p.format..), method =  "wilcox.test"),label.x = c(1), label.y = c(73))  +
  scale_y_break(breaks=c(22,72), scales = "free")+
  facet_wrap(~Malignant_clusters, scales = "free_y") + theme_classic() + theme(text=element_text(size=15))+NoLegend()
dev.copy2pdf(file="MALIGNANT/MES_Cell_percentages_RvsNR.pdf", width = 4, height = 5,useDingbats=FALSE,family="sans")



#cell_percentage_no <-  cell_percentage %>% subset(Malignant_clusters %not_in% c("Antigen_presenting", "Mesenchymal_like"))
cell_percentage_no <-  cell_percentage %>% subset(Malignant_clusters %not_in% c("Patient_specific_A", "Patient_specific_B"))
ggboxplot(cell_percentage_no, x = "Timepoint", y = "percentage",
         fill  = "Timepoint", palette =c("#377eb8", "#4daf4a"),
          shape = "Timepoint",
          add = "jitter")+
  stat_compare_means(aes(group = Timepoint, label = paste0("p = ", ..p.format..), method =  "wilcox.test"),label.x = c(1))  +
  facet_wrap(~Malignant_clusters, nrow = 1, scales = "free_y") + theme_classic() + theme(text=element_text(size=10))
dev.copy2pdf(file="MALIGNANT/Cell_percentages_TIMEPOINT.pdf", width = 16, height = 3,useDingbats=FALSE,family="sans")


#response + timepoint
cell_percentage_no <-  cell_percentage %>% subset(Malignant_clusters %not_in% c("Antigen_Presentation", "Mesenchymal_like"))
cell_percentage_no <-  cell_percentage_no %>% subset(Malignant_clusters %not_in% c("Patient_specific_A", "Patient_specific_B"))

cell_percentage_no$var <- paste(cell_percentage_no$Timepoint,cell_percentage_no$Response, sep = "_" )
table(cell_percentage_no$var)
my_comp <- list(c("OT_NR", "OT_R"), c("BT_NR", "BT_R"), c("BT_R", "OT_R"), c("BT_NR", "OT_NR"))
ggboxplot(cell_percentage_no, x = "var", y = "percentage",
          fill = "Timepoint",
          shape = "Response",
          palette = c("#377eb8","#4daf4a"),
          add = "jitter")+
  stat_compare_means(comparisons =  my_comp, label = "p.format", method =  "wilcox.test") +
  facet_wrap(~Malignant_clusters,nrow = 1, scales = "free_y") + theme_classic() + RotatedAxis() +  theme(text=element_text(size=12),axis.text.x = element_text(angle = 45,hjust = 1))
dev.copy2pdf(file="MALIGNANT/Cell_percentages_RvsNR_Timepoint.pdf",  width = 16, height = 5,useDingbats=FALSE,family="sans")

cell_percentage_no <-  cell_percentage %>% subset(Malignant_clusters %in% c("Antigen_Presentation"))
cell_percentage_no$var <- paste(cell_percentage_no$Timepoint,cell_percentage_no$Response, sep = "_" )
table(cell_percentage_no$var)
my_comp <- list(c("OT_NR", "OT_R"), c("BT_NR", "BT_R"), c("BT_R", "OT_R"), c("BT_NR", "OT_NR"))
ggboxplot(cell_percentage_no, x = "var", y = "percentage",
          fill = "Timepoint",
          shape = "Response",
          palette = c("#377eb8","#4daf4a"),
          add = "jitter")+
  stat_compare_means(comparisons =  my_comp, label = "p.format", method =  "wilcox.test") +
  facet_wrap(~Malignant_clusters,nrow = 1, scales = "free_y") + theme_classic() + RotatedAxis() +  theme(text=element_text(size=15),axis.text.x = element_text(angle = 45,hjust = 1)) + NoLegend()
dev.copy2pdf(file="MALIGNANT/Antigen_Cell_percentages_RvsNR_Timepoint.pdf", width = 4, height = 5,useDingbats=FALSE,family="sans")


cell_percentage_no <-  cell_percentage %>% subset(Malignant_clusters %in% c("Mesenchymal_like"))
cell_percentage_no$var <- paste(cell_percentage_no$Timepoint, cell_percentage_no$Response, sep = "_" )
table(cell_percentage_no$var)
my_comp <- list(c("OT_NR", "OT_R"), c("BT_NR", "BT_R"), c("BT_R", "OT_R"), c("BT_NR", "OT_NR"))
ggboxplot(cell_percentage_no, x = "var", y = "percentage",
          fill = "Timepoint",
          shape = "Response",
          palette = c("#377eb8","#4daf4a"),
          add = "jitter")+
  stat_compare_means(comparisons =  my_comp, label = "p.format", method =  "wilcox.test") +
  scale_y_break(breaks=c(22,72))+
  facet_wrap(~Malignant_clusters,nrow = 1) + theme_classic() + RotatedAxis() +  theme(text=element_text(size=15),axis.text.x = element_text(angle = 45,hjust = 1))+ NoLegend()
dev.copy2pdf(file="MALIGNANT/Mes_Cell_percentages_RvsNR_Timepoint.pdf", width = 4, height = 5,useDingbats=FALSE,family="sans")


pdf("MALIGNANT/Cell_percentages_TILs_bin.pdf", width = 19, height = 8)
subset_tils <- cell_percentage %>% subset(TILs_bin %in% c("Present", "Absent"))
ggboxplot(subset_tils, x = "TILs_bin", y = "percentage",
          fill  = "TILs_bin",
          shape = "TILs_bin",
          add = "jitter")+
  stat_compare_means(aes(group = TILs_bin, label = paste0("p = ", ..p.format..), method =  "wilcox.test"),label.x = c(2)) +
  facet_wrap(~Malignant_clusters, nrow = 2, scale = "free") + theme_classic() + theme(text=element_text(size=18),axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

ggboxplot(subset_tils, x = "Response", y = "percentage",
          color = "Response",
          shape = "TILs_bin",
          add = "jitter")+
  stat_compare_means(aes(group = TILs_bin, label = paste0("p = ", ..p.format..), method =  "wilcox.test"),label.x = c(1), label.y = 75) +
  facet_wrap(~Malignant_clusters, nrow = 2) + theme_classic() + theme(text=element_text(size=18),axis.text.x = element_text(angle = 45, hjust = 1))

pdf("MALIGNANT/Cell_percentages_TILs.pdf", width = 15, height = 8)
subset_tils <- cell_percentage %>% subset(TILs %in% c("Brisk", "NonBrisk", "Absent", "Necrosis"))
ggboxplot(subset_tils, x = "TILs", y = "percentage",
          fill = "TILs",
          shape = "TILs",
          add = "jitter")+
  stat_compare_means(aes(group = TILs, label = paste0("p = ", ..p.format..), method =  "kruskal.test"), label.x = c(1.5), label.y = 75) +
  facet_wrap(~Malignant_clusters, nrow = 2,  scales = "free_y") + theme_classic() + theme(text=element_text(size=15),axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()



pdf("MALIGNANT/Cell_percentages_Tissue.pdf", width = 13, height = 5)
ggboxplot(cell_percentage, x = "Tissue", y = "percentage",
          color = "Tissue",
          #shape = "Tissue",
          add = "jitter")+
  stat_compare_means(aes(group = Tissue, label = paste0("p = ", ..p.format..), method =  "kruskal.test"), label.x = c(1.5), label.y = 75) +
  facet_wrap(~Malignant_clusters, nrow = 2) + theme_classic() + theme(text=element_text(size=15),axis.text.x = element_text(angle = 45,hjust = 1))
dev.off()

pdf("MALIGNANT/Cell_percentages_Mutation.pdf", width = 13, height = 6)
ggboxplot(cell_percentage, x = "Mutation_type", y = "percentage",
          fill = "Mutation_type",
          #shape = "Mutation_type",
          add = "jitter")+
  stat_compare_means(aes(group = Mutation_type, label = paste0("p = ", ..p.format..), method =  "kruskal.test"), label.x = c(1.5)) +
  facet_wrap(~Malignant_clusters, nrow = 1, scale = "free") + theme_classic() + theme(text=element_text(size=12),axis.text.x = element_text(angle = 45,hjust = 1))
dev.off()



ggplot(cell_percentage, aes(orig.ident, y=percentage, fill=Malignant_clusters))+geom_bar(stat="identity", position = "stack")  +  facet_grid(~GC_number, scale="free", drop = TRUE) + theme_classic()+RotatedAxis()
RES <- subset(cell_percentage, Response == "R")
pdf("MALIGNANT/Cell_percentages_bar_chart_R.pdf", width = 15, height = 12)
ggplot(RES, aes(orig.ident, y=percentage, fill=Malignant_clusters))+geom_bar(stat="identity", position = "stack")  +  facet_grid(~GC_number, scale="free", drop = TRUE) + theme_classic()+RotatedAxis()
dev.off()
NRES <- subset(cell_percentage, Response == "NR")
pdf("MALIGNANT/Cell_percentages_bar_chart_NR.pdf", width = 15, height = 12)
ggplot(NRES, aes(orig.ident, y=percentage, fill=Malignant_clusters))+geom_bar(stat="identity", position = "stack")  +  facet_grid(~GC_number+Timepoint+Response, scale="free", drop = TRUE) + theme_classic()+RotatedAxis()
dev.off()
pdf("MALIGNANT/Cell_percentages_bar_chart_ALL_SAMPLES.pdf", width = 24, height = 12)
ggplot(cell_percentage, aes(orig.ident, y=percentage, fill=Malignant_clusters))+geom_bar(stat="identity", position = "stack")  +  facet_grid(~GC_number+Timepoint+Response, scale="free", drop = TRUE) + theme_classic()+RotatedAxis()
dev.off()



################paired

table(cell_percentage$Timepoint,  cell_percentage$GC_number)
table(cell_percentage$GC_number,  cell_percentage$orig.ident)
### subset the samples which do not have pairs
`%not_in%` <- purrr::negate(`%in%`)
subset_for_matched <- cell_percentage %>% subset(GC_number %not_in% c("22",
                                                                      "25",
                                                                      "28",
                                                                      "29",
                                                                      "32",
                                                                      "37", "38", "39"))
subset_for_matched <- subset_for_matched %>% subset(orig.ident %not_in% c("scrCMA050" )) # patient 16 select randomly sample from multiple biospsies from the same patient
subset_for_matched <- subset_for_matched %>%
  mutate(GC_number = ifelse(orig.ident == "scrCMA040", "16_bis", GC_number))
subset_for_matched <- subset_for_matched %>%
  mutate(GC_number = ifelse(orig.ident == "scrCMA048", "16_bis", GC_number))

########## MES
cell_percentage_no <-  subset_for_matched %>% subset(Malignant_clusters %in% c( "Mesenchymal_like"))
cell_percentage_no$var <- paste(cell_percentage_no$Timepoint,cell_percentage_no$Response, sep = "_" )
table(cell_percentage_no$var)
my_comp <- list(c("BT_NR", "OT_NR"), c("BT_R", "OT_R"))
ggboxplot(cell_percentage_no, x = "var", y = "percentage",
          fill = "Timepoint",
          shape = "Response",
          add = "jitter",
          palette = c("#377eb8","#4daf4a"),
          id = "GC_number")+ ylim(0,80)+
  stat_compare_means(comparisons =  my_comp, method =  "wilcox.test", paired= T, exact = F, label.y = c(73,73)) +
  geom_line(aes(group = GC_number), color = "gray", alpha = 0.4) + 
    scale_y_break(c(22,72))+
  facet_wrap(~Malignant_clusters,nrow = 1, scales = "free_y") + theme_classic() + RotatedAxis() +  theme(text=element_text(size=12),axis.text.x = element_text(angle = 45,hjust = 1))
dev.copy2pdf(file="MALIGNANT/Mes_Cell_percentages_RvsNR_Timepoint_PAIRED.pdf", width = 5, height = 5,useDingbats=FALSE,family="sans")

########### APC
cell_percentage_no <-  subset_for_matched %>% subset(Malignant_clusters %in% c( "Antigen_presentation"))

cell_percentage_no$var <- paste(cell_percentage_no$Timepoint,cell_percentage_no$Response, sep = "_" )
table(cell_percentage_no$var)
my_comp <- list(c("BT_NR", "OT_NR"), c("BT_R", "OT_R"))
ggboxplot(cell_percentage_no, x = "var", y = "percentage",
          fill = "Timepoint",
          shape = "Response",
          add = "jitter",
          palette = c("#377eb8","#4daf4a"),
          id = "GC_number")+ ylim(0,45)+
  stat_compare_means(comparisons =  my_comp, method =  "wilcox.test", paired= T, exact = F, label.y = c(42,42)) +
  geom_line(aes(group = GC_number), color = "gray", alpha = 0.4) + 
  scale_y_break(c(40,41))+
  facet_wrap(~Malignant_clusters,nrow = 1, scales = "free_y") + theme_classic() + RotatedAxis() +  theme(text=element_text(size=12),axis.text.x = element_text(angle = 45,hjust = 1))
dev.copy2pdf(file="MALIGNANT/APC_Cell_percentages_RvsNR_Timepoint_PAIRED.pdf", width = 5, height = 5,useDingbats=FALSE,family="sans")



#### only aPD1
cell_percentage_aPD1 <-  cell_percentage %>% subset(Malignant_clusters %not_in% c("Patient_specific_A", "Patient_specific_B"))
cell_percentage_aPD1 <-  cell_percentage_aPD1 %>% subset(Treatment %not_in% c("NNivolumabIpilimumab"))
cell_percentage_aPD1$var <- paste(cell_percentage_aPD1$Timepoint,cell_percentage_aPD1$Response, sep = "_" )
table(cell_percentage_aPD1$var)
my_comp <- list(c("OT_NR", "OT_R"), c("BT_NR", "BT_R"), c("BT_R", "OT_R"), c("BT_NR", "OT_NR"))
ggboxplot(cell_percentage_aPD1, x = "var", y = "percentage",
          fill = "Timepoint",
          shape = "Response",
          palette = c("#377eb8","#4daf4a"),
          add = "jitter")+
  stat_compare_means(comparisons =  my_comp, label = "p.format", method =  "wilcox.test") +
  facet_wrap(~Malignant_clusters,nrow = 1, scales = "free_y") + theme_classic() + RotatedAxis() +  theme(text=element_text(size=10),axis.text.x = element_text(angle = 45,hjust = 1))
dev.copy2pdf(file="MALIGNANT/Cell_percentages_RvsNR_Timepoint_onlyPD1.pdf", width = 15, height = 5,useDingbats=FALSE,family="sans")






