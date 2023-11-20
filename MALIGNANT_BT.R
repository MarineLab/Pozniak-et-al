library(nichenetr) #used for converting mouse to human gene symbols
GC_all_malignant<- readRDS("/MALIGNANT/GC_ALL_MALIGNANT.rds")

GC_all_malignant@meta.data$"Malignant_clusters" <- plyr::revalue(as.character(GC_all_malignant$seurat_clusters),
                                                                 c( "0" = "Melanocytic",
                                                                    "1" = "Mitochondrial(low_quality)",
                                                                    "2" = "Melanocytic",
                                                                    "3" = "Antigen_presentation",
                                                                    "4" = "Interferon_alpha_beta_response",
                                                                    "5" ="Stress (p53 response)",
                                                                    "6" ="Neural_Crest_like",
                                                                    "7" ="Stress (hypoxia response)",
                                                                    "8" ="Mitotic",
                                                                    "9" ="Patient_specific_A",
                                                                    "10" ="Mesenchymal_like",
                                                                    "11" ="Patient_specific_B"
                                                                 ))
DimPlot(GC_all_malignant, split.by = "BT/OT")
table(GC_all_malignant$`GC number`)

GC_malignant_subset_BT <- subset(GC_all_malignant, subset = `BT/OT` == "BT")
GC_malignant_subset_BT
table(GC_malignant_subset_BT$`GC number`)
table(GC_malignant_subset_BT$`orig.ident`)

pdf("/mean_cnv_per_mt.pdf", width = 9, height = 5)
VlnPlot(GC_malignant_subset_BT, features = c("mean_cnv"), group.by = "Malignant_clusters",pt.size = 0)+  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))
VlnPlot(GC_malignant_subset_BT, features = c("percent.mt"), group.by = "Malignant_clusters",pt.size = 0)+  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))
dev.off()

Idents(GC_malignant_subset_BT) <- "seurat_clusters"
GC_malignant_subset_BT_markers <- FindAllMarkers(GC_malignant_subset_BT, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.3) # 0.3, 0.4 before
write.table(GC_malignant_subset_BT_markers, "/GC_malignant_subset_BT_markers_seurat_clusters.txt", sep='\t', quote = FALSE, col.names = T, row.names = F)

Idents(GC_malignant_subset_BT) <- "Malignant_clusters"

####################################### Cluster Identity #########################
GC_malignant_subset_BT_markers <- FindAllMarkers(GC_malignant_subset_BT, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.3) # 0.3, 0.4 before
write.table(GC_malignant_subset_BT_markers, "/GC_malignant_subset_BT_markers.txt", sep='\t', quote = FALSE, col.names = T, row.names = F)
top10 <- GC_malignant_subset_BT_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top5 <- GC_malignant_subset_BT_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)

pdf("/Heatmap_GC_malignant_subset_BT.pdf", width = 15, height = 15)
DoHeatmap(GC_malignant_subset_BT, features = top10$gene) + NoLegend()
dev.off()
pdf("/Heatmap_GC_malignant_subset_BT_multibar_candidate_genes.pdf", width = 15, height = 15)
DoMultiBarHeatmap(GC_malignant_subset_BT, 
                               features = genes, 
                               cells = NULL, 
                               group.by = "Malignant_clusters", 
                               additional.group.by = "orig.ident") 
dev.off()

Idents(GC_malignant_subset_BT) <- "Malignant_clusters"

cluster.averages <- AverageExpression(GC_malignant_subset_BT, return.seurat = TRUE)
mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)

my_levels <- c("Antigen_presentation",  "Interferon_alpha_beta_response","Melanocytic","Mesenchymal_like", "Mitochondrial(low_quality)","Mitotic","Neural_Crest_like", "Patient_specific_A", "Patient_specific_B", "Stress (hypoxia response)", "Stress (p53 response)")
levels(cluster.averages) <- my_levels

Idents(GC_all_malignant) <- "Malignant_clusters"
levels(GC_all_malignant)

# [1] "Neural_Crest_like"              "Interferon_alpha_beta_response" "Melanocytic"                    "Mesenchymal_like"               "Mitochondrial(low_quality)"    
#[6] "Antigen_presentation"           "Mitotic"                        "Stress (p53 response)"          "Stress (hypoxia response)"      "Patient_specific_B"            
#[11] "Patient_specific_A"     
palette.malignant <- c("#77a600","#08519c","#ba5e45","#836ab1","grey", "#ff7f00" ,"#b30024", "#fcbd1f",  "#f28ac6", "#20b2aa","#e02887")

#pdf("/UMAP_Malignant_clusters.pdf", width = 6, height = 4)
DimPlot(GC_malignant_subset_BT, group.by = c("Malignant_clusters"), cols = palette.malignant, pt.size = 1)+ NoAxes()
#dev.off()
dev.copy2pdf(file="/UMAP_Malignant_clusters.pdf", width = 7, height = 4,useDingbats=FALSE,family="sans")

pdf("/UMAP_seurat_clusters.pdf", width = 4, height = 4)
DimPlot(GC_malignant_subset_BT)+ NoAxes()
dev.off()
pdf("/UMAP_orig.ident.pdf", width = 8, height = 7)
DimPlot(GC_malignant_subset_BT, group.by = c("orig.ident"), pt.size = 1)+ NoAxes()
dev.off()
pdf("/UMAP_GC.pdf", width = 8, height = 7)
DimPlot(GC_malignant_subset_BT, group.by = c("GC number"), pt.size = 1)+ NoAxes()
dev.off()
pdf("/UMAP_mut.pdf", width = 8, height = 7)
DimPlot(GC_malignant_subset_BT, group.by = c("Mutation"), pt.size = 1)+ NoAxes()
dev.off()
GC_malignant_subset_BT@meta.data$"TILs" <- plyr::revalue(as.character(GC_malignant_subset_BT$orig.ident),
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
###################### percentages  out of malignant cells ##################
###########################################################################

GC_malignant_subset_BT@meta.data$"TILs_bin" <- plyr::revalue(as.character(GC_malignant_subset_BT$TILs),
                                                       c("Absent" = "Absent",
                                                         "NonBrisk" = "Present",
                                                         "Brisk" = "Present"
                                                       ))
GC_malignant_subset_BT@meta.data$"Mutation" <- plyr::revalue(as.character(GC_malignant_subset_BT$`Mut type`),
                                                       c("N/A" = "NA"
                                                       ))
GC_malignant_subset_BT@meta.data$"Tissue" <- plyr::revalue(as.character(GC_malignant_subset_BT$`Tissue`),
                                                     c("Lymph node" = "LN"
                                                     ))

library(readxl)
Meta_data <- read_excel("/Meta_data_GEX_RESPONSE_2.0-25-08-2021.xlsx")
Meta_data <- as.data.frame(Meta_data)
GC_malignant_subset_BT@meta.data$row_names <- rownames(GC_malignant_subset_BT@meta.data)
dim(GC_malignant_subset_BT@meta.data)
GC_malignant_subset_BT@meta.data<-GC_malignant_subset_BT@meta.data %>% left_join(Meta_data, by="orig.ident")
row.names(GC_malignant_subset_BT@meta.data) <- GC_malignant_subset_BT@meta.data$row_names
table(GC_malignant_subset_BT$orig.ident)

TEST <- GC_malignant_subset_BT@meta.data
cell_num <- TEST %>%
  mutate(sample_id = as.factor(paste(orig.ident,`BT/OT`, `GC number`, `TILs_bin`, `TILs`, `Tissue`, `Mutation`,`Response`,  sep="_"))) %>%
  mutate(Malignant_clusters = as.factor(Malignant_clusters)) %>%
  group_by(sample_id, Malignant_clusters, .drop=FALSE) %>%
  dplyr::summarise(n=n()) %>%
  tidyr::separate(sample_id, c("orig.ident", "Timepoint", "GC_number", "TILs_bin", "TILs", "Tissue", "Mutation_type", "Response"))

total_cells<- TEST %>%
  group_by(orig.ident) %>%
  dplyr::summarise(total = n())
total_cells
cell_percentage<- left_join(cell_num, total_cells) %>%
  mutate(percentage = n/total*100)
cell_percentage
  cell_percentage <- cell_percentage %>%
  mutate(Response = recode(Response , "0" = "NR", "1" = "R"))

pdf("/Cell_percentage_response.pdf", width = 14, height = 7)
ggboxplot(cell_percentage, x = "Response", y = "percentage",
          color = "Response",
          shape = "Response",
          add = "jitter")+
  stat_compare_means(aes(group = Response, label = paste0("p = ", ..p.format..), method =  "wilcox.test"), label.x = c(1), label.y = 80) +
  facet_wrap(~Malignant_clusters, ncol = 5) + theme_classic() + theme(text=element_text(size=20),axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

pdf("/Cell_percentages_TILs_bin.pdf", width = 17, height = 5)
subset_tils <- cell_percentage %>% subset(TILs_bin %in% c("Present", "Absent"))
ggboxplot(subset_tils, x = "TILs_bin", y = "percentage",
          fill = "TILs_bin",
          shape = "TILs_bin",
          add = "jitter")+
  stat_compare_means(aes(group = TILs_bin, label = paste0("p = ", ..p.format..), method =  "wilcox.test"), label.x = c(1), label.y = 80) +
  facet_wrap(~Malignant_clusters, nrow = 1) + theme_classic() + theme(text=element_text(size=20),axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

pdf("/Cell_percentages_TILs.pdf", width = 17, height = 5)
my_comp <- list(c("Brisk", "Absent"), c("Brisk", "NonBrisk"), c("NonBrisk", "Absent"))

subset_tils <- cell_percentage %>% subset(TILs %in% c("Brisk", "NonBrisk", "Absent"))
subset_tils <-  subset_tils %>% subset(Malignant_clusters %not_in% c("Patient_specific_A", "Patient_specific_B"))

ggboxplot(subset_tils, x = "TILs", y = "percentage",
          fill  = "TILs",
          shape = "TILs",
          add = "jitter")+
  stat_compare_means(comparisons =  my_comp, label = "p.format", method =  "wilcox.test") + NoLegend()+
  facet_wrap(~Malignant_clusters, nrow = 1) + theme_classic() + theme(text=element_text(size=9),axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


pdf("/Cell_percentages_Tissue.pdf", width = 14, height = 7)
ggboxplot(cell_percentage, x = "Tissue", y = "percentage",
          color = "Tissue",
          #shape = "Tissue",
          add = "jitter")+
  stat_compare_means(aes(group = Tissue, label = paste0("p = ", ..p.format..), method =  "kruskal.test"), label.x = c(1), label.y = 90) +
  facet_wrap(~Malignant_clusters, ncol = 5) + theme_classic() + theme(text=element_text(size=20),axis.text.x = element_text(angle = 45,hjust = 1))
dev.off()

pdf("/Cell_percentages_Mutation.pdf", width = 7, height = 6)
cell_percentage_no_NA <- cell_percentage %>% subset(Mutation_type !="NA")
cell_percentage_no_NA <- cell_percentage_no_NA %>% subset(Mutation_type !="WT")

cell_percentage_no_NA <-  cell_percentage_no_NA %>% subset(Malignant_clusters %not_in% c("Patient_specific_A", "Patient_specific_B"))

ggboxplot(cell_percentage_no_NA, x = "Mutation_type", y = "percentage",
          fill = "Mutation_type",
          #shape = "Mutation_type",
          add = "jitter")+
  stat_compare_means(aes(group = Mutation_type, label = paste0("p = ", ..p.format..), method =  "kruskal.test"), label.x = c(1.5), label.y = 80) +
  facet_wrap(~Malignant_clusters, ncol = 3) + theme_classic() + theme(text=element_text(size=14),axis.text.x = element_text(angle = 45,hjust = 1))
dev.off()

pdf("/Stackbar_per_patient.pdf", width = 10, height = 5)
#cell_percentage_malig <- cell_percentage %>% dplyr::select(Malignant_clusters, percentage,orig.ident,Response, Timepoint) %>% tidyr::spread(Malignant_clusters, percentage)
custom_colors <- list()
colourCount = 11
getPalette = colorRampPalette(brewer.pal(11, "Dark2"))
Malignant_clusters <- getPalette(colourCount)
ggplot(cell_percentage, aes(orig.ident, y=percentage, fill=Malignant_clusters))+geom_bar(stat="identity", position = "stack") +  theme(text=element_text(size=40))+  scale_fill_manual(values = getPalette(colourCount)) + theme_classic()+RotatedAxis()
dev.off()

pdf("/Stackbar_per_Patient.pdf", width = 10, height = 5)
#cell_percentage_malig <- cell_percentage %>% dplyr::select(Malignant_clusters, percentage,orig.ident,Response, Timepoint) %>% tidyr::spread(Malignant_clusters, percentage)
ggplot(cell_percentage, aes(orig.ident, y=percentage, fill=Malignant_clusters))+geom_bar(stat="identity", position = "stack") +  facet_grid(~GC_number, scales = "free_x", space = "free_x", drop = TRUE) +  theme(text=element_text(size=40)) + scale_fill_manual(values = getPalette(colourCount)) + theme_classic()+RotatedAxis()
dev.off()

pdf("/Stackbar_per_patient_mutation.pdf", width = 10, height = 5)
#cell_percentage_malig <- cell_percentage %>% dplyr::select(Malignant_clusters, percentage,orig.ident,Response, Timepoint) %>% tidyr::spread(Malignant_clusters, percentage)
ggplot(cell_percentage, aes(orig.ident, y=percentage, fill=Malignant_clusters))+geom_bar(stat="identity", position = "stack", preserve = "single") +  facet_grid(~Mutation_type, scales = "free_x", space = "free_x", drop = TRUE) +  theme(text=element_text(size=40)) + theme_classic()+RotatedAxis()
dev.off()
pdf("/n_cells_per_patient.pdf", width = 10, height = 5)
ggplot(GC_malignant_subset_BT@meta.data, aes(orig.ident))+geom_bar(stat="count") + theme(text=element_text(size=40)) + theme_classic()+RotatedAxis()
dev.off()

levels(GC_malignant_subset_BT$Malignant_clusters) <- c("Antigen_presentation", "Interferon_alpha_beta","Melanocytic_OXPHOS", "Mesenchymal", "Mitochondrial","Mitotic", "Neural_like", "Patient_specific_A", "Patient_specific_B", "Stress(hypoxia)", "Trans_reg")
Idents(GC_malignant_subset_BT) <- "Malignant_clusters"


##########################        HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION        ##########################
gmtFile <- paste("/HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.gmt")
geneSets <- getGmt(gmtFile)
geneSets
cells_rankings <- AUCell_buildRankings(GC_malignant_subset_BT@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION<-getAUC(cells_AUC)
HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION<-t(HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION)
GC_malignant_subset_BT@meta.data<-cbind(GC_malignant_subset_BT@meta.data, HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION)
pdf("/HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.pdf", width = 7.08, height = 5.8)
FeaturePlot(GC_malignant_subset_BT, features = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", label = T) 
dev.off()
pdf("/HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_Vln.pdf", width = 7.08, height = 5)
VlnPlot(GC_malignant_subset_BT, features = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", pt.size = 0, group.by = "Malignant_clusters") + theme(plot.margin=unit(c(1,1,1.5,1.5),"cm"), legend.position = 'none')
dev.off()
FeatureScatter(GC_malignant_subset_BT, "PRRX1", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")

##########################        Rambow_MESlike        ##########################
genes<-c("ADM","ANGPTL4","BCAT1","CCL2","CDH2","CYSLTR2","DLC1","DLX1","EDNRA","ERRFI1","FABP4","FGF1","GPC3","IGFBP6","LMO4","LOX","LOXL2","MGP","NDNF","NES","NR2F1","PDGFRB","PRDX1","PTGER4","RGS16","SH2B3","SLIT2","SPRY2","TGM2","TMSB4X","UNC5B","VEGFA","VSNL1","IL13RA2","FOSL2","AXL","PLXDC1","CDH13","COL1A2","COL3A1","RGS5","VCAN","IGFBP5","DDAH1","TGFBI","SOX4","BGN","COL1A1","TNC")
geneSets <- GeneSet(genes, setName="Rambow_MESlike")
geneSets
cells_rankings <- AUCell_buildRankings(GC_malignant_subset_BT@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Rambow_MESlike<-getAUC(cells_AUC)
Rambow_MESlike<-t(Rambow_MESlike)
GC_malignant_subset_BT@meta.data<-cbind(GC_malignant_subset_BT@meta.data, Rambow_MESlike)
pdf("/Rambow_MESlike.pdf", width = 7.08, height = 5.8)
FeaturePlot(GC_malignant_subset_BT, features = "Rambow_MESlike", label = T)
dev.off()
pdf("/Rambow_MESlike_Vln.pdf", width = 5, height = 4)
VlnPlot(GC_malignant_subset_BT, features = "Rambow_MESlike", pt.size = 0, group.by = "Malignant_clusters") + theme(plot.margin=unit(c(1,1,1.5,1.5),"cm"), legend.position = 'none')
dev.off()

##########################        Rambow_SMC        ##########################
genes<-c("C7ORF53","ATP6AP1L","UBXN10","CD36","DLX5","ARMC7","PAX3","IP6K3","TRIM67","SLC7A8","SLC25A48","RNF121","KIAA1161","FRAT2")
geneSets <- GeneSet(genes, setName="Rambow_SMC")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_malignant_subset_BT@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Rambow_SMC<-getAUC(cells_AUC)
Rambow_SMC<-t(Rambow_SMC)
GC_malignant_subset_BT@meta.data<-cbind(GC_malignant_subset_BT@meta.data, Rambow_SMC)
pdf("/Rambow_SMC.pdf", width = 7.08, height = 5.8)
FeaturePlot(GC_malignant_subset_BT, features = "Rambow_SMC", label = T)
dev.off()
pdf("/Rambow_SMC_Vln.pdf", width = 5, height = 4)
VlnPlot(GC_malignant_subset_BT, features = "Rambow_SMC", pt.size = 0, group.by = "Malignant_clusters")  + theme(plot.margin=unit(c(1,1,1.5,1.5),"cm"), legend.position = 'none')
dev.off()

##########################        Rambow_Hyperdiff        ##########################
genes<-c("TYRP1",  "SLC24A5",  "PMEL",  "KIT",  "GPR143",  "APOE",  "MLANA",  "TRPM1",  "MLPH",  "TYR",  "SNAI2",  "SLC45A2",  "EDNRB",  "RAB27A",  "FABP7")
geneSets <- GeneSet(genes, setName="Rambow_Hyperdiff")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_malignant_subset_BT@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Rambow_Hyperdiff<-getAUC(cells_AUC)
Rambow_Hyperdiff<-t(Rambow_Hyperdiff)
GC_malignant_subset_BT@meta.data<-cbind(GC_malignant_subset_BT@meta.data, Rambow_Hyperdiff)
pdf("/Rambow_Hyperdiff.pdf", width = 7.08, height = 5.8)
FeaturePlot(GC_malignant_subset_BT, features = "Rambow_Hyperdiff", label = T)
dev.off()
pdf("/Rambow_Hyperdiff_Vln.pdf", width = 5, height = 4)
VlnPlot(GC_malignant_subset_BT, features = "Rambow_Hyperdiff", pt.size = 0, group.by = "Malignant_clusters")  + theme(plot.margin=unit(c(1,1,1.5,1.5),"cm"), legend.position = 'none')
dev.off()
##########################          Rambow_NCSC       ##########################
gene5 <- c("A2M",  "ADAMTS4",  "C6ORF103",  "ANXA1",  "AQP1",  "ATP1A2",  "ATP1B2",  "CADM1",  "CNN3",  "COL1A1",  "COL4A1",  "GDNF",  "GFRA1",  "GFRA2",  "GFRA3",  "IGF1",  "IL1RAP",  "ITGA1",  "ITGA6",  "L1CAM",  "LAMC1",  "MATN2",  "MEF2C",  "MPZ",  "NGFR",  "NLGN3",  "NRXN1",  "PDGFB",  "PLAT",  "PRIMA1",  "RSPO3",  "RXRG",  "S100A4",  "SEMA3B",  "SLC22A17",  "SLITRK6",  "SOX10",  "SYT11",  "TFAP2B",  "THBS2",  "TMEM176B",  "VCAN")
geneSets <- GeneSet(gene5, setName="Rambow_NCSC")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_malignant_subset_BT@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Rambow_NCSC<-getAUC(cells_AUC)
Rambow_NCSC<-t(Rambow_NCSC)
GC_malignant_subset_BT@meta.data<-cbind(GC_malignant_subset_BT@meta.data, Rambow_NCSC)
pdf("/Rambow_NCSC.pdf", width = 7.08, height = 5.8)
FeaturePlot(GC_malignant_subset_BT, features = "Rambow_NCSC", label = T)
dev.off()
pdf("/Rambow_NCSC_vln.pdf", width = 5, height = 4)
VlnPlot(GC_malignant_subset_BT, features = c("Rambow_NCSC"), pt.size = 0, group.by = "Malignant_clusters") + theme(plot.margin=unit(c(1,1,1.5,1.5),"cm"), legend.position = 'none')
dev.off()

##########################        Tirosh_AXL         ##########################
genes<-c("ANGPTL4", "FSTL3", "GPC1", "TMSB10", "SH3BGRL3", "PLAUR", "NGFR", "SEC14L2", "FOSL1", "SERPINE1", "IGFBP3", "TNFRSF12A", "GBE1", "AXL", "PHLDA2", "MAP1B", "GEM", "SLC22A4", "TYMP", "TREM1", "RIN1", "S100A4", "COL6A2", "FAM46A", "CITED1", "S100A10", "UCN2", "SPHK1", "TRIML2", "S100A6", "TMEM45A", "CDKN1A", "UBE2C", "ERO1L", "SLC16A6", "CHI3L1", "FN1", "S100A16", "CRIP1", "SLC25A37", "LCN2", "ENO2", "PFKFB4", "SLC16A3", "DBNDD2", "LOXL2", "CFB", "CADM1", "LTBP3", "CD109", "AIM2", "TCN1", "STRA6", "C9orf89", "DDR1", "TBC1D8", "METTL7B", "GADD45A", "UPP1", "SPATA13", "GLRX", "PPFIBP1", "PMAIP1", "COL6A1", "JMJD6", "CIB1", "HPCAL1", "MT2A", "ZCCHC6", "IL8", "TRIM47", "SESN2", "PVRL2", "DRAP1", "MTHFD2", "SDC4", "NNMT", "PPL", "TIMP1", "RHOC", "GNB2", "PDXK", "CTNNA1", "CD52", "SLC2A1", "BACH1", "ARHGEF2", "UBE2J1", "CD82", "ZYX", "P4HA2", "PEA15", "GLRX2", "HAPLN3", "RAB36", "SOD2", "ESYT2", "IL18BP", "FGFRL1", "PLEC")
geneSets <- GeneSet(genes, setName="Tirosh_AXL")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_malignant_subset_BT@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Tirosh_AXL<-getAUC(cells_AUC)
Tirosh_AXL<-t(Tirosh_AXL)
GC_malignant_subset_BT@meta.data<-cbind(GC_malignant_subset_BT@meta.data, Tirosh_AXL)
pdf("/Tirosh_AXL.pdf", width = 5, height = 4) 
FeaturePlot(GC_malignant_subset_BT, features = "Tirosh_AXL", label = T) 
dev.off()
pdf("/Tirosh_AXL_vln.pdf", width = 5, height =4) 
VlnPlot(GC_malignant_subset_BT, features = "Tirosh_AXL", pt.size = 0, group.by = "Malignant_clusters") + theme(plot.margin=unit(c(1,1,1.5,1.5),"cm"), legend.position = 'none')
dev.off()

##########################        Tirosh_MITF        ##########################
genes<-c("MITF", "TYR", "PMEL", "PLP1", "GPR143", "MLANA", "STX7", "IRF4", "ERBB3", "CDH1", "GPNMB", "IGSF11", "SLC24A5", "SLC45A2", "RAP2B", "ASAH1", "MYO10", "GRN", "DOCK10", "ACSL3", "SORT1", "QPCT", "S100B", "MYC", "LZTS1", "GYG2", "SDCBP", "LOXL4", "ETV5", "C1orf85", "HMCN1", "OSTM1", "ALDH7A1", "FOSB", "RAB38", "ELOVL2", "MLPH", "PLK2", "CHL1", "RDH11", "LINC00473", "RELL1", "C21orf91", "SCAMP3", "SGK3", "ABCB5", "SLC7A5", "SIRPA", "WDR91", "PIGS", "CYP27A1", "TM7SF3", "PTPRZ1", "CNDP2", "CTSK", "BNC2", "TOB1", "CELF2", "ROPN1", "TMEM98", "CTSA", "LIMA1", "CD99", "IGSF8", "FDFT1", "CPNE3", "SLC35B4", "EIF3E", "TNFRSF14", "VAT1", "HPS5", "CDK2", "CAPN3", "SUSD5", "ADSL", "PIGY", "PON2", "SLC19A1", "KLF6", "MAGED1", "ERGIC3", "PIR", "SLC25A5", "JUN", "ARPC1B", "SLC19A2", "AKR7A2", "HPGD", "TBC1D7", "TFAP2A", "PTPLAD1", "SNCA", "GNPTAB", "DNAJA4", "APOE", "MTMR2", "ATP6V1B2", "C16orf62", "EXOSC4", "STAM")
geneSets <- GeneSet(genes, setName="Tirosh_MITF")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_malignant_subset_BT@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Tirosh_MITF<-getAUC(cells_AUC)
Tirosh_MITF<-t(Tirosh_MITF)
GC_malignant_subset_BT@meta.data<-cbind(GC_malignant_subset_BT@meta.data, Tirosh_MITF)
pdf("/Tirosh_MITF.pdf", width = 5, height = 4)
FeaturePlot(GC_malignant_subset_BT, features = "Tirosh_MITF", label = T) 
dev.off()
pdf("/Tirosh_MITF_vln.pdf", width = 5, height = 4)
VlnPlot(GC_malignant_subset_BT, features = "Tirosh_MITF", pt.size = 0, group.by = "Malignant_clusters") + theme(plot.margin=unit(c(1,1,1.5,1.5),"cm"), legend.position = 'none')
dev.off()


##########################        Hoek_INV        ##########################
genes<-c("ADAM12", "AMOTL2", "AXL", "BIRC3", "CDH13", "CDK14", "COL13A1", "CRIM1", "CRISPLD2", "CYR61", "DPYD", "EFEMP1", "EGFR", "F2RL1", "FGF2", "FLNB", "FOXD1", "FST", "FZD2", "HEG1", "HS3ST3A1", "ITGA2", "ITGA3", "KCNMA1", "LOXL2", "MYOF", "NRP1", "NTM", "NUAK1", "OSMR", "PDGFC", "PODXL", "S100A2", "SLC22A4", "SLIT2", "SYNJ2", "TCF4", "THBS1", "TLE4", "TNFRSF11B", "TPBG", "TPM1", "TRAM2", "WNT5A", "ZEB1")
geneSets <- GeneSet(genes, setName="Hoek_INV")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_malignant_subset_BT@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Hoek_INV<-getAUC(cells_AUC)
Hoek_INV<-t(Hoek_INV)
GC_malignant_subset_BT@meta.data<-cbind(GC_malignant_subset_BT@meta.data, Hoek_INV)
pdf("/Hoek_INV.pdf", width = 4, height = 5)
FeaturePlot(GC_malignant_subset_BT, features = "Hoek_INV",label = T)  + theme(legend.position = 'none')
dev.off()
pdf("/HOEK_INVASIVE_vln.pdf", width = 5, height = 4)
VlnPlot(GC_malignant_subset_BT, features = "Hoek_INV", pt.size = 0, group.by = "Malignant_clusters") + theme(plot.margin=unit(c(1,1,1.5,1.5),"cm"), legend.position = 'none')
dev.off()

##########################        Hoek_PRO         ##########################
genes<-c("ACP5", "ADCY2", "APOE", "ASAH1", "BIRC7", "C21orf91", "CAPN3", "CDH1", "CDK2", "CDK5R1", "CEACAM1", "DAPK1", "DCT", "FAM174B", "GALNT3", "GNPTAB", "GPM6B", "GPR143", "GPRC5B", "GYG2", "HPS4", "INPP4B", "IRF4", "IVNS1ABP", "KAZN", "MBP", "MICAL1", "MITF", "MLANA", "MYO1D", "NR4A3", "OCA2", "PHACTR1", "PIR", "PLXNC1", "PMEL", "RAB27A", "RAB38", "RGS20", "RHOQ", "RRAGD", "SEMA6A", "SIRPA", "SLC45A2", "ST3GAL6", "STX7", "TNFRSF14", "TRPM1", "TYR", "TYRP1", "WDR91", "ZFYVE16")
geneSets <- GeneSet(genes, setName="Hoek_PRO")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_malignant_subset_BT@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Hoek_PRO<-getAUC(cells_AUC)
Hoek_PRO<-t(Hoek_PRO)
GC_malignant_subset_BT@meta.data<-cbind(GC_malignant_subset_BT@meta.data, Hoek_PRO)
pdf("/Hoek_PRO.pdf", width = 7.08, height = 5.8)
FeaturePlot(GC_malignant_subset_BT, features = "Hoek_PRO",label = T)
dev.off()
pdf("/Hoek_PRO_vln.pdf", width = 5, height = 4)
VlnPlot(GC_malignant_subset_BT, features = "Hoek_PRO", pt.size = 0, group.by = "Malignant_clusters") + theme(plot.margin=unit(c(1,1,1.5,1.5),"cm"), legend.position = 'none')
dev.off()




##########################          Tsoi_undiff       ##########################
gene5 <- c("AJUBA", "TOR4A", "MARCH4", "ZDHHC2", "ZNF467", "ZNF185", "ZIC2", "VASN", "UCP2", "GALNT6", "TNFAIP2", "TNFSF18", "TMEM40", "TMEM200A", "TMEM184A", "TBL1X", "TRERF1", "TOX", "TBC1D2", "SFN", "SAMD12", "SAMD11", "SOX9", "SLC8A1", "SLC38A4", "SLC16A14", "SCN5A", "SCNN1A", "SH3RF2", "SERPINB7", "SLPI", "SECTM1", "RUNX2", "ARHGAP29", "REN", "PAWR", "PSG9", "PSG5", "PSG4", "PBX1", "PLAGL1", "PHLDB2", "PLEKHA6", "PDGFC", "PLAU", "PKP2", "PLAC8", "PADI3", "PITX1", "NUAK1", "NTNG1", "NMT2", "MYEOV", "MICAL2", "MGST1", "MECOM", "LYPD6B", "LAMA5", "KISS1", "KRT86", "KRT81", "KRT80", "KRT8", "KRT7", "KRT18", "JUP", "IL7R", "IL4R", "IRS1", "IGFN1", "HES7", "GDA", "GLIS2", "GATA2", "GPRC5C", "GPRC5A", "FMNL1", "FOXA1", "FLNC", "FERMT1", "FAT4", "FAM196B", "ELFN2", "EGFR", "DSE", "DMBT1", "DIO2", "DOCK2", "CYP2S1", "CRIM1", "CDK15", "CORO6", "COLEC10", "CCDC88C", "CCDC69", "F3", "F2RL1", "CLU", "CDYL2", "CITED2", "CARD11", "CPA4", "CREB3L1", "CNN1", "CALB2", "CDH4", "BTBD11", "BDNF", "BASP1", "BNC1", "ATP8B1", "ABCG2", "ARMC4", "ANKRD1", "AR", "AMIGO2", "ADAMTSL1", "ACSL5")
geneSets <- GeneSet(gene5, setName="Tsoi_undiff")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_malignant_subset_BT@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Tsoi_undiff<-getAUC(cells_AUC)
Tsoi_undiff<-t(Tsoi_undiff)
GC_malignant_subset_BT@meta.data<-cbind(GC_malignant_subset_BT@meta.data, Tsoi_undiff)
pdf("/Tsoi_undiff.pdf", width = 7.08, height = 5.8)
FeaturePlot(GC_malignant_subset_BT, features = "Tsoi_undiff", label = T) 
dev.off()
pdf("/Tsoi_undiff_vln.pdf", width = 5, height = 4)
VlnPlot(GC_malignant_subset_BT, features = c("Tsoi_undiff"), pt.size = 0, group.by = "Malignant_clusters") + theme(plot.margin=unit(c(1,1,1.5,1.5),"cm"), legend.position = 'none')
dev.off()


##########################          Tsoi NC_like       ##########################
gene5 <- c("PXYLP1",  "CXCL8",  "CEMIP",  "TCAF2",  "ZNF469",  "WNT5A",  "TMEM47",  "TMEM171",  "TGFBI",  "TGFA",  "TFAP2C",  "TSPAN13",  "SQRDL",  "SULF1",  "ST8SIA5",  "SOX2",  "SLC24A3",  "SLITRK6",  "SHISA2",  "SH3PXD2A",  "SERTAD4",  "STK32B",  "SEMA3B",  "SFRP1",  "S100A6",  "RAMP1",  "PMEPA1",  "PCSK5",  "PHLDA2",  "PLA2G7",  "OPRD1",  "NTM",  "NRXN3",  "NES",  "MUC5B",  "MAP1LC3A",  "LRRC15",  "KIAA1755",  "ITGB8",  "IER3",  "HHEX",  "GDNF",  "GLI2",  "FOXC2",  "FLT1",  "FAT3",  "FEZ1",  "FAM135B",  "EHF",  "EML1",  "DRD2",  "DEPDC7",  "CYB5R2",  "CSRP2",  "CCL2",  "CADM3",  "CADM1",  "CD96",  "CTSS",  "CHST2",  "CHST1",  "CACNA2D3",  "BST1",  "ABCA6",  "ANGPTL4",  "AIM2")
geneSets <- GeneSet(gene5, setName="Tsoi_NClike")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_malignant_subset_BT@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Tsoi_NClike<-getAUC(cells_AUC)
Tsoi_NClike<-t(Tsoi_NClike)
GC_malignant_subset_BT@meta.data<-cbind(GC_malignant_subset_BT@meta.data, Tsoi_NClike)
pdf("/Tsoi_NClike.pdf", width = 7.08, height = 5.8)
FeaturePlot(GC_malignant_subset_BT, features = "Tsoi_NClike", label = T)
dev.off()
pdf("/Tsoi_NClike_vln.pdf", width = 5, height = 4)
VlnPlot(GC_malignant_subset_BT, features = c("Tsoi_NClike"), pt.size = 0, group.by = "Malignant_clusters") + theme(plot.margin=unit(c(1,1,1.5,1.5),"cm"), legend.position = 'none')
dev.off()


##########################          Tsoi Transitory       ##########################
gene5 <- c("XYLT1", "TSPAN7", "SOD3", "SCRG1", "SORL1", "SEMA3E", "SELENBP1", "RNASE1", "RAPGEF4", "PCDH7", "PRSS33", "PCSK6", "PLBD1", "NELL1", "NPR1", "MCAM", "MMP15", "MAMDC2", "LSAMP", "LRRTM4", "GDF11", "FXYD3", "EBF3", "COL11A2", "COL9A1", "CX3CL1", "BCHE", "ANO4", "ALDH1A1")
geneSets <- GeneSet(gene5, setName="Tsoi_Transit")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_malignant_subset_BT@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Tsoi_Transit<-getAUC(cells_AUC)
Tsoi_Transit<-t(Tsoi_Transit)
GC_malignant_subset_BT@meta.data<-cbind(GC_malignant_subset_BT@meta.data, Tsoi_Transit)
pdf("/Tsoi_Transit.pdf", width = 7.08, height = 5.8)
FeaturePlot(GC_malignant_subset_BT, features = "Tsoi_Transit", label = T)
dev.off()
pdf("/Tsoi_Transit_vln.pdf", width = 5, height = 4)
VlnPlot(GC_malignant_subset_BT, features = c("Tsoi_Transit"), pt.size = 0, group.by = "Malignant_clusters") + theme(plot.margin=unit(c(1,1,1.5,1.5),"cm"), legend.position = 'none')
dev.off()


##########################          Tsoi Melanocytic       ##########################
gene5 <- c("CCDC171", "CFAP61", "ZDHHC11B", "VEPH1", "TNFRSF14", "TDRD3", "TPPP", "TRIM63", "TRPM1", "TTC39A", "TSPAN10", "SLC7A8", "SLC16A6", "SLAMF7", "SEMA6A", "RUNX3", "RNF144B", "RNLS", "RGS12", "PYCARD", "PRUNE2", "PRKCB", "PRDM7", "KCNAB2", "OCA2", "NR4A3", "NAV2", "MYO1D", "MAPK4", "MAT1A", "MLANA", "LXN", "KCP", "IL16", "IL12RB2", "HSD17B14", "HMOX1", "H2AFJ", "GOLGA7B", "QPCT", "GFOD1", "GPR143", "FYB", "FAM83H", "FAM174B", "EPHA5", "ENTHD1", "DNAJA4", "DENND2D", "C2orf88", "CCL18", "CEACAM1", "CAPG", "CDH3", "CDH1", "ATP6V0D2", "ABCD1", "ABCB5", "APOLD1", "ANKRD30B", "ADCY2", "ADAM23")
geneSets <- GeneSet(gene5, setName="Tsoi_Melanocytic")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_malignant_subset_BT@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Tsoi_Melanocytic<-getAUC(cells_AUC)
Tsoi_Melanocytic<-t(Tsoi_Melanocytic)
GC_malignant_subset_BT@meta.data<-cbind(GC_malignant_subset_BT@meta.data, Tsoi_Melanocytic)
pdf("/Tsoi_Melanocytic.pdf", width = 7.08, height = 5.8)
FeaturePlot(GC_malignant_subset_BT, features = "Tsoi_Melanocytic", label = T)
dev.off()
pdf("/Tsoi_Melanocytic_vln.pdf", width = 5, height =4)
VlnPlot(GC_malignant_subset_BT, features = c("Tsoi_Melanocytic"), pt.size = 0, group.by = "Malignant_clusters") + theme(plot.margin=unit(c(1,1,1.5,1.5),"cm"), legend.position = 'none')
dev.off()

##########################        BARON_Melanocytic        ##########################
genes<-c("ACTB",  "ACTR1B",  "ALDH7A1",  "ALDH9A1",  "ANP32A",  "ANXA1",  "ANXA2",  "ANXA4",  "ANXA5",  "AP2M1",  "AQP3",  "ARPC3",  "ARPC5",  "ARRDC3",  "ATP5B",  "ATP5C1",  "ATP5D",  "ATP5F1",  "ATP5G1",  "ATP5O",  "ATP6V0E1",  "BLOC1S6",  "BTG3",  "C15orf48",  "C21orf33",  "CALM1",  "CALML3",  "CAPG",  "CAPZB",  "CCDC80",  "CCT3",  "CCT8",  "CD81",  "CDC42SE1",  "CFL1",  "CIRBP",  "CKAP2",  "CKB",  "CLU",  "CNBP",  "COBLL1",  "COPB2",  "COPE",  "COPS4",  "COX4I2",  "COX6B1",  "COX7A2",  "CSRP2",  "CSTB",  "CTSL",  "DAD1",  "DAP",  "DAPL1",  "DBI",  "DCT",  "DDC",  "DDOST",  "DENR",  "DYNC1LI2",  "EEF1G",  "EEF2",  "EFEMP2",  "EHD1",  "EIF4A1",  "EIF5A2",  "ELOVL1",  "ENO1",  "ETF1",  "FHL1",  "FNDC5",  "FSTL1",  "FXYD6-FXYD2",  "GAPDHS",  "GDI1",  "GLRX",  "GM2A",  "GPR143",  "GPX3",  "GPX4",  "GSTP1",  "GSTT1",  "GYG1",  "H2AFZ",  "HIGD1A",  "HINT1",  "HM13",  "HMGA1",  "HMGB2",  "HMGN2",  "HNRNPA0",  "HNRNPD",  "HSPE1",  "IFI44",  "KARS",  "KMT5A",  "KRT8",  "LDHB",  "LMNB1",  "LMNB2",  "LOXL2",  "LSM7",  "MARCKSL1",  "MDH1",  "MDH2",  "MIF",  "MTPN",  "MYH7B",  "MYL12A",  "MYL6",  "NDUFA1",  "NDUFA8",  "NPM1",  "PABPC1",  "PAH",  "PAICS",  "PALM",  "PBDC1",  "PCBD1",  "PFN1",  "PFN2",  "PKM",  "PLPP1",  "PMEL",  "POMP",  "POR",  "PPIA",  "PPIB",  "PPP2CB",  "PRDX1",  "PRDX2",  "PRDX6",  "PSMA2",  "PSMA3",  "PSMB6",  "PSMC3",  "PSMC5",  "PSME2",  "PSMG1",  "PTGR1",  "QDPR",  "RAD21",  "RAN",  "RAP1B",  "REEP5",  "RHOC",  "RPN1",  "RTN1",  "S100A10",  "S100A4",  "SEC13",  "SEC61B",  "SELENOF",  "SLC25A5",  "SLC39A10",  "SOD1",  "SPON2",  "SSR4",  "SUB1",  "SYNGR3",  "TCP1",  "TFG",  "TIMP2",  "TKTL2",  "TMA7",  "TMED2",  "TMEM163",  "TTR",  "TUBB4B",  "TYRP1",  "UBE2L3",  "UFM1",  "UQCRC2",  "VIM",  "YWHAB",  "YWHAE",  "YWHAQ")
geneSets <- GeneSet(genes, setName="BARON_Melanocytic")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_malignant_subset_BT@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
BARON_Melanocytic<-getAUC(cells_AUC)
BARON_Melanocytic<-t(BARON_Melanocytic)
GC_malignant_subset_BT@meta.data<-cbind(GC_malignant_subset_BT@meta.data, BARON_Melanocytic)
pdf("/BARON_Melanocytic.pdf", width = 7.08, height = 5.8)
FeaturePlot(GC_malignant_subset_BT, features = "BARON_Melanocytic",label = T)
dev.off()
pdf("/BARON_Mature_vln.pdf", width = 5, height = 4)
VlnPlot(GC_malignant_subset_BT, features = "BARON_Melanocytic", pt.size = 0, group.by = "Malignant_clusters") + theme(plot.margin=unit(c(1,1,1.5,1.5),"cm"), legend.position = 'none')
dev.off()

##########################        BARON_STRESS        ##########################
genes<-c("AHSA1",  "ATF3",  "B2M",  "BRD2",  "BTG1",  "BTG2",  "CCER1",  "CCNG1",  "CHCHD2",  "CIRBP",  "CYR61",  "DNAJB1",  "DUSP1",  "DUSP2",  "DUSP5",  "EGR1",  "EGR2",  "FOS",  "FOSB",  "FOSL1",  "GADD45B",  "HES1",  "HES4",  "HIST2H2AB",  "HSP90AA1",  "HSPA1L",  "HSPA4",  "HSPA5",  "HSPA8",  "HSPH1",  "ID2",  "ID3",  "IER2",  "IER5",  "JDP2",  "JUN",  "JUNB",  "JUND",  "MCL1",  "MIDN",  "MKNK2",  "NFKBIA",  "NR4A1",  "NR4A3",  "PDCD4",  "PLK3",  "PNRC2",  "PPP1R15A",  "RHBDD1",  "RSRP1",  "SGK1",  "SKIL",  "SLC38A2",  "SOX4",  "SQSTM1",  "STIP1",  "TOB1",  "UBB",  "UBC",  "WSB1",  "ZFAND2B")
geneSets <- GeneSet(genes, setName="BARON_STRESS")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_malignant_subset_BT@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
BARON_STRESS<-getAUC(cells_AUC)
BARON_STRESS<-t(BARON_STRESS)
GC_malignant_subset_BT@meta.data<-cbind(GC_malignant_subset_BT@meta.data, BARON_STRESS)
pdf("/BARON_STRESS.pdf", width = 7.08, height = 5.8)
FeaturePlot(GC_malignant_subset_BT, features = "BARON_STRESS",label = T)
dev.off()
pdf("/BARON_STRESS_vln.pdf", width = 5, height = 4)
VlnPlot(GC_malignant_subset_BT, features = "BARON_STRESS", pt.size = 0, group.by = "Malignant_clusters") + theme(plot.margin=unit(c(1,1,1.5,1.5),"cm"), legend.position = 'none')
dev.off()

##########################        BARON_Neural_like        ##########################
genes<-c("ABCF1",  "ABI2",  "ACSL4",  "ACTN1",  "ACVRL1",  "ADCYAP1",  "ADD3",  "AK2",  "AKAP8L",  "ALCAM",  "AP1B1",  "APLP2",  "ARL6IP1",  "ATAD2B",  "ATF4",  "ATP1B1",  "ATP1B3",  "ATP5A1",  "ATP5E",  "ATP5G3",  "BCKDK",  "BLCAP",  "BTF3",  "C6orf62",  "CALM3",  "CAPN2",  "CAPZA1",  "CAV1",  "CCND1",  "CCNI",  "CDH1",  "CDK15",  "CDKN1B",  "CEP170",  "CFH",  "CHUK",  "CIART",  "CLIC4",  "CNTFR",  "COL9A2",  "COX6C",  "CSNK2A1",  "CSNK2B",  "CSTB",  "CTSH",  "CTSL",  "CTSZ",  "CXXC5",  "CYC1",  "CYCS",  "DCBLD1",  "DDX21",  "DUSP6",  "EDNRA",  "EEF1A2",  "EEF1B2",  "EEF1D",  "EIF2S2",  "EIF2S3",  "EIF3G",  "EIF3I",  "EIF4B",  "EIF4EBP2",  "EIF5",  "EIF5A",  "EMILIN1",  "EMP3",  "ENPP1",  "FABP3",  "FAM107B",  "FAM126A",  "FAU",  "FGFR4",  "FNBP1L",  "GADD45A",  "GLO1",  "GLUL",  "GNB1",  "GNG2",  "GNG3",  "GNPDA2",  "GPI",  "H1FX",  "H2AFY2",  "HDAC1",  "HDLBP",  "HNRNPAB",  "HNRNPR",  "HPS1",  "HSPA9",  "HSPB1",  "IFI27",  "IFI27L1",  "IFI27L2",  "ILF2",  "ITGA5",  "ITM2B",  "JAK1",  "JAM3",  "KPNB1",  "LGMN",  "LYPLA2",  "METAP1",  "MITF",  "MLPH",  "MOV10",  "MPEG1",  "MT2A",  "MTMR6",  "MYCBP",  "MYCN",  "MYO5B",  "NACA",  "NCL",  "NDRG1",  "NDUFA4",  "NDUFS5",  "NME1",  "NMRK2",  "NPEPL1",  "NSA2",  "NTNG2",  "NTRK3",  "NUDC",  "OAZ1",  "OLA1",  "PA2G4",  "PABPC1",  "PABPN1",  "PDCD6IP",  "PDLIM1",  "PEF1",  "PFDN5",  "PGAM1",  "PGK1",  "PIK3CA",  "PLA2G7",  "PLEC",  "PNO1",  "PNOC",  "PPP1R1B",  "PRMT1",  "PSAP",  "PSMB1",  "PTP4A2",  "PTPRE",  "PVALB",  "QKI",  "RAB5A",  "RABL6",  "RACK1",  "RDH8",  "RGL1",  "RGS12",  "RHPN2",  "RNASET2",  "ROMO1",  "RPL22L1",  "RPLP2",  "RPS26",  "RSL24D1",  "RTN4",  "RUNX3",  "SASH1",  "SDC4",  "SDHA",  "SDPR",  "SEH1L",  "40057",  "SERBP1",  "SF3A3",  "SGPL1",  "SH3BGRL3",  "SH3GLB1",  "SHFM1",  "SLC20A1",  "SLC25A22",  "SLC2A11",  "SMAP1",  "SNCG",  "SNRPB",  "SNRPD1",  "SOX10",  "SOX2",  "SPIRE1",  "SPRY4",  "SSB",  "ST8SIA6",  "STAT2",  "SUMO1",  "SYT11",  "TBC1D10A",  "TCN2",  "TFAP2A",  "TFAP2C",  "TFAP2E",  "TFRC",  "TMC5",  "TMEM30A",  "TMEM50A",  "TMEM59L",  "TP53I11",  "TPT1",  "TXNIP",  "U2AF2",  "UQCRC1",  "VCP",  "VPS26A",  "WBP2NL",  "YWHAB",  "YWHAG",  "ZBTB20")
geneSets <- GeneSet(genes, setName="BARON_NClike")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_malignant_subset_BT@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
BARON_NClike<-getAUC(cells_AUC)
BARON_NClike<-t(BARON_NClike)
GC_malignant_subset_BT@meta.data<-cbind(GC_malignant_subset_BT@meta.data, BARON_NClike)
pdf("/BARON_NClike.pdf", width = 7.08, height = 5.8)
FeaturePlot(GC_malignant_subset_BT, features = "BARON_NClike",label = T)
dev.off()
pdf("/BARON_NClike.pdf", width = 5, height = 4)
VlnPlot(GC_malignant_subset_BT, features = "BARON_NClike", pt.size = 0, group.by = "Malignant_clusters") + theme(plot.margin=unit(c(1,1,1.5,1.5),"cm"), legend.position = 'none')
dev.off()



##########################        Verfaillie_PRO        ##########################
genes<-c("GDPD5",  "SLC22A23",  "RP11-3L8.3",  "PDZRN3",  "NR4A1",  "PROS1",  "SS18L1",  "AP1S2",  "TSPAN10",  "BHLHE41",  "CCDC171",  "IL17D",  "PDK4",  "C2orf88",  "PRR5",  "CDK2",  "CTSH",  "MREG",  "FAM124A",  "RP11-558F24.4",  "ASAH1",  "HEY2",  "SDC3",  "RAP1GAP",  "CAPG",  "RAB6B",  "SORT1",  "MXI1",  "PCSK9",  "ST3GAL4",  "APOC1",  "ADAMTS17",  "EPHA5",  "HEY1",  "SLC12A7",  "SH3TC1",  "PAG1",  "RAB38",  "RP11-1055B8.2",  "SCIN",  "LINC00340",  "GAB2",  "CES3",  "DUSP15",  "CYP27A1",  "ATRNL1",  "FAXC",  "HES6",  "KLF15",  "DGCR5",  "AC005786.5",  "INPP5F",  "SGK1",  "NPAS1",  "CYGB",  "PIR",  "CERS1",  "AATK",  "SESN3",  "CDK5R1",  "SCUBE2",  "SLC17A9",  "SLC16A10",  "TMTC2",  "KIAA1598",  "RASSF2",  "MAST1",  "RP11-80F22.9",  "TBC1D7",  "PRKCH",  "PDE3B",  "C19orf71",  "SLAIN1",  "GJA3",  "FAM53B",  "RAB11FIP4",  "BAI1",  "C11orf96",  "SLC27A3",  "FGF13",  "BRSK2",  "EGLN3",  "GNAL",  "CABLES1",  "GPR137B",  "CXADR",  "SHC2",  "ST6GALNAC1",  "FBXL16",  "Z83851.1",  "ASRGL1",  "TNFRSF19",  "CEACAM1",  "SORL1",  "ANKRD6",  "ISG20",  "RIMS4",  "MYOM2",  "LAD1",  "ADRBK2",  "LZTS1",  "RNF125",  "TMEM255A",  "TRPM8",  "FAM20A",  "LONRF3",  "PPM1H",  "SPTBN2",  "TEX41",  "ST6GAL1",  "POU3F3",  "ADAM23",  "ANO4",  "MFSD12",  "RP11-390P2.4",  "ST6GALNAC2",  "RP11-137H2.6",  "OLFM2",  "TMCC2",  "GREB1",  "TTC39A",  "FAM213A",  "KCNS1",  "TNFRSF14",  "STXBP6",  "ITGA7",  "ALDH1A1",  "ZNF704",  "BAMBI",  "PGBD5",  "PRKCZ",  "IL6R",  "PLCL1",  "EGR3",  "ITPKB",  "NAT16",  "LRRC4",  "STOX2",  "OGDHL",  "PIK3AP1",  "PNLIPRP3",  "CNTN3",  "BAAT",  "COL25A1",  "CELF2",  "RASIP1",  "TMEM229B",  "PLEKHG1",  "PKNOX2",  "KRTAP19-1",  "SLC7A4",  "SLC24A4",  "ASB4",  "ST3GAL6-AS1",  "EFR3B",  "NKAIN1",  "CASKIN1",  "LDLRAD4",  "WNK2",  "TKTL1",  "RAB17",  "KNDC1",  "TESK2",  "CHN2",  "SLC7A8",  "FGD4",  "CRTAC1",  "PPARGC1A",  "RP11-527H14.2",  "FAM134B",  "RASEF",  "ACAN",  "CHST6",  "TFCP2L1",  "HMCN1",  "RP3-395M20.8",  "PNMAL1",  "RP3-527G5.1",  "PIP5K1B",  "IL12RB2",  "TENM1",  "RPS6KA2",  "CECR2",  "VGF",  "BCAN",  "ADCY1",  "RAB3C",  "CLDN14",  "RP11-481A20.11",  "MERTK",  "LINC00937",  "LINC00504",  "CCL18",  "RXRG",  "PHACTR1",  "CARD14",  "QPCT", 
         "GRASP",  "DLL3",  "PKLR",  "LRGUK",  "C1orf51",  "TEX15",  "B4GALNT3",  "KIAA1211",  "PELI2",  "POU3F2",  "PRKCB",  "HCG20",  "RP11-98L5.2",  "DISC1FP1",  "RP11-317M11.1",  "KBTBD11",  "LAMC3",  "RP11-557H15.4",  "MPZ",  "AC011294.3",  "RENBP",  "PRUNE2",  "LAMA1",  "B3GAT1",  "TINCR",  "NDN",  "TTYH2",  "CTB-151G24.1",  "TC2N",  "SLC16A6",  "RP11-143A12.3",  "LINC00518",  "MFI2",  "PRODH",  "GOLGA7B",  "POU3F1",  "CDH3",  "GPM6A",  "NR4A3",  "MAPT",  "LARGE",  "LRP2",  "LPL",  "TMPRSS5",  "RLBP1",  "GNG7",  "ACP5",  "RP11-93B14.5",  "FAM167B",  "EDNRB",  "TMC6",  "LINC00426",  "CA8",  "MYO16",  "ST3GAL6",  "CITED1",  "RASGEF1A",  "AC009784.3",  "DAPK1",  "KIT",  "GYG2",  "PLEKHB1",  "AP000479.1",  "TMPRSS13",  "NUP210",  "ALDH1A2",  "ZNF536",  "FAM19A5",  "RGS1",  "NKX2-5",  "OPLAH",  "TSPAN7",  "KDR",  "PLEKHH1",  "SBK1",  "FAM155B",  "ITGA9",  "BEST1",  "KCNAB2",  "RP11-2E17.1",  "MAGEB2",  "RP13-735L24.1",  "TUBB8P7",  "CPVL",  "DENND1C",  "MOB3B",  "ST8SIA6",  "WIPF3",  "PRSS33",  "CNTN1",  "CCDC64",  "APOD",  "NRG3",  "RP3-395M20.7",  "AC009499.1",  "PPP1R14C",  "CAPN3",  "SOX6",  "MBP",  "ISM1",  "MITF",  "CACNA1H",  "RP4-718J7.4",  "MGAT4A",  "CPB2-AS1",  "LRRC4B",  "EXTL1",  "SOX8",  "CRYAB",  "ROPN1B",  "SYT3",  "SGCA",  "NAT8L",  "RP11-509E16.1",  "PMP2",  "NMRK2",  "AC004988.1",  "GSTT1",  "LINC00589",  "RP11-189B4.6",  "TYRP1",  "FAM189A2",  "ALDH3B2",  "RP11-161M6.2",  "RAB33A",  "ZDHHC11B",  "ROBO2",  "SLC35F1",  "RRAGD",  "LINGO1",  "RP11-290F20.3",  "HSPB8",  "RP11-669N7.2",  "RP3-332B22.1",  "CACNA1D",  "SAMD5",  "GFPT2",  "ITGAX",  "CPN1",  "RP11-726G1.1",  "MSI1",  "TNRC6C-AS1",  "LGI3",  "MLIP",  "TUBB4A",  "COBL",  "LHFPL3-AS1",  "LINC00488",  "RP11-557H15.2",  "GALNT3",  "S100B",  "AC002511.1",  "SORBS1",  "IGF1",  "ADCY2",  "DLGAP1",  "SFTPC",  "HILS1",  "PLXNC1",  "MYH14",  "AC145110.1",  "PLA1A",  "GLB1L2",  "CHL1",  "GAPDHS",  "CTNNA2",  "RP11-599J14.2",  "SHE",  "CTD-2207A17.1",  "RP11-1055B8.3",  "RP11-104E19.1",  "AC096559.1",  "FXYD3",  "HPGD",  "MMP8",  "IRX6",  "RP11-347E10.1",  "GJB1",  "GPR143",  "PLP1",  "CDH19",  "IL16",  "MPPED2",  "RP11-252C15.1",  "CDH1",  "C10orf90",  "APOE",  "ITIH5",  "OCA2",  
         "MAF", "FAM69C", "RP4-529N6.1",  "SOX10",  "FRMD4B",  "ERBB3",  "SEMA6A",  "MCF2L",  "MAPK4",  "FAM174B",  "SLC38A8",  "MYO1D",  "LINC00520",  "NSG1",  "IGSF11",  "C20orf26",  "HRK",  "PAEP",  "SLC45A2",  "TRIM51",  "ATP10A",  "ROPN1",  "COL9A3",  "SCML4",  "RP11-429E11.2",  "GAS7",  "CA14",  "LINC00698",  "PTPRZ1",  "DCT",  "MID2",  "PMEL",  "SORCS1",  "ABCB5",  "PRDM7",  "GPM6B",  "LCP2",  "ENTHD1",  "MLANA",  "TYR",  "NKAIN4",  "IRF4",  "SGCD",  "SLC24A5",  "TRIM63",  "BIRC7",  "TRPM1")
geneSets <- GeneSet(genes, setName="Verfaillie_PRO")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_malignant_subset_BT@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Verfaillie_PRO<-getAUC(cells_AUC)
Verfaillie_PRO<-t(Verfaillie_PRO)
GC_malignant_subset_BT@meta.data<-cbind(GC_malignant_subset_BT@meta.data, Verfaillie_PRO)
pdf("/Verfaillie_PRO.pdf", width = 7.08, height = 5.8)
FeaturePlot(GC_malignant_subset_BT, features = "Verfaillie_PRO",label = T)
dev.off()
pdf("/Verfaillie_PRO_vln.pdf", width = 5, height = 4)
VlnPlot(GC_malignant_subset_BT, features = "Verfaillie_PRO", pt.size = 0, group.by = "Malignant_clusters") + theme(plot.margin=unit(c(1,1,1.5,1.5),"cm"), legend.position = 'none')
dev.off()


##########################        Verfaillie_INV        ##########################
genes<-c("MXRA5",  "GABRE",  "DDX3Y",  "VCAM1",  "POSTN",  "CD248",  "ABCC3",  "ANPEP",  "TMEM200A",  "ADAMTS12",  "F2RL2",  "IL1B",  "RPS4Y1",  "EREG",  "COL1A1",  "SERPINB7",  "ADAMTSL1",  "KDM5D",  "TXLNG2P",  "GPR128",  "ALPK2",  "SLC14A1",  "VIPR1",  "CDH13",  "PTGFR",  "PTX3",  "CPA4",  "SERPINE1",  "CCL2",  "CCDC68",  "NETO1",  "EDN1",  "ELTD1",  "PLAU",  "RP11-48O20.4",  "CDH6",  "EIF1AY",  "RGS4",  "COL8A1",  "GPR39",  "FENDRR",  "ESM1",  "CTA-392C11.1",  "DPEP3",  "SSTR1",  "PRKY",  "FPR3",  "SMEK3P",  "TNFRSF9",  "PXDN",  "FPR1",  "MAGEA4", 
         "FAM196B",  "COL1A2",  "ANXA10",  "MKX",  "WNT7B",  "ARMC4",  "FGF5",  "MUC13",  "UTY",  "GALNT5",  "GDF5",  "LRRC15",  "KAL1",  "SPANXD",  "LRRC17",  "SMOC1",  "NTNG1",  "BDNF",  "IFI27",  "ZFY",  "PAPPA",  "LOXL2",  "IL1A",  "USP9Y",  "TSPYL5",  "SPANXC",  "CCBE1",  "LYPD6B",  "PAPPA2",  "ADAMTS6",  "C7orf69",  "CPA3",  "IGFN1",  "C3",  "G0S2",  "TRHDE",  "ITGA11",  "CDA",  "CTD-2171N6.1",  "ABCC9",  "APBB1IP",  "IL6",  "FAM43B",  "ITGBL1",  "CREB3L1",  "TTTY15",  "FOXF1",  "FOXR2",  "PSG5",  "VGLL3",  "FOXG1",  "TNFSF18",  "RAB27B",  "IDO1",  
         "FAM180A",  "PARM1",  "RAB3B",  "KRT81",  "MEOX2",  "GREM1",  "ECSCR",  "GALNT6",  "IL32",  "CD163L1",  "TRHDE-AS1",  "C6orf141",  "MMP19",  "THBD",  "LDB2",  "PBX1",  "GPR68",  "BDKRB2",  "TENM2",  "CFI",  "EDIL3",  "LINC00707",  "MIR137HG",  "GBP1",  "F2RL1",  "PDE1C",  "RTN1",  "AXL",  "ABI3BP",  "TRIML2",  "NEXN",  "CA9",  "TCF4",  "NID2",  "MYCT1",  "HS3ST3A1",  "CLMP",  "SLC15A3",  "GSTM1",  "SH2D4A",  "COL13A1",  "GLIS3",  "BNC1",  "POU2F2",  "SPATA18",  "COL5A1",  "AOX1",  "S1PR1",  "SCN9A",  "MT1E",  "EFEMP1",  "GLIPR1",  "GBP4",  "RP11-371I1.2",  
         "OAS2",  "NMNAT2",  "RP11-224O19.2",  "RSAD2",  "CHMP4C",  "FAM155A",  "PITX1",  "CFH",  "GABRQ",  "TRIM58",  "PRDM8",  "LOXL1",  "APOL3",  "ARHGDIB",  "JPH2",  "KRT80",  "EMILIN1",  "TNFRSF11B",  "HTR1F",  "RRAD",  "IL18R1",  "IL7R",  "AMIGO2",  "ATP8B1",  "FBN2",  "CXCL2",  "NRG1",  "ISG15",  "HSPB6",  "INHBA",  "RASGRF2",  "NTN4",  "GUCY1B3",  "IL4I1",  "COL12A1",  "MT1A",  "IL11",  "TMEM158",  "STK33",  "ARNTL2",  "PKP2",  "LAMC2",  "TOR4A",  "SOD3",  "TMEM119",  "SRGN",  "IGFBP6",  "SCG2",  "IL1RL1",  "TOX2",  "NLRP3",  "BIRC3",  "HSPB7",  "FBN1",  "SULF1", 
         "KCNMA1",  "PROSER2-AS1",  "SSC5D",  "DENND2A",  "PDCD1LG2",  "S1PR3",  "COL6A3",  "LPAR1",  "VIT",  "COL3A1",  "PTPRN",  "DNAH2",  "DKK3",  "FGF1",  "XAF1",  "HECW2",  "RP11-166D19.1",  "CFB",  "FLNC",  "ACE",  "CYP1B1",  "GALNT13",  "SLC6A10P",  "COL5A2",  "THBS1",  "SERPINB2",  "CLIC2",  "MPP4",  "KLHL4",  "CA12",  "RARRES3",  "VEGFC",  "COL4A6",  "WNT5A",  "HHIPL2",  "OASL",  "CLDN11",  "CSF1R",  "MYPN",  "S100A4",  "TGM2",  "SLIT2",  "TNFAIP2",  "KRT18",  "CCL5",  "IFI44L",  "MAGEA10",  "MX2",  "NRP1",  "CCDC80",  "MAGEA11",  "GABRA2",  "HRH1",  "IFITM1",  "NR2F1",  
         "TRPC4",  "TGFB2",  "CDK15",  "ZNF385D",  "KRT8",  "RP5-1011O1.3",  "OXTR",  "NGEF",  "LOX",  "FMOD",  "PPARG",  "RP11-513O17.2",  "SRPX2",  "IRAK2",  "MGST1",  "CRISPLD2",  "AR",  "EPHB2",  "LAMA3",  "CDH12",  "TRIM22",  "RP11-400K9.4",  "PTPLAD2",  "ARHGAP22",  "CHRDL1",  "SLFN11",  "ADAM19",  "COLEC12",  "STAC",  "MYL9",  "NTSR1",  "TSLP",  "PDGFRB",  "STC2",  "PCDHGA7",  "GPR1",  "SLC16A9",  "EGFR",  "SPOCD1",  "FRMPD4",  "AKR1C3",  "PRTFDC1",  "CLU",  "FAM101A",  "PLAGL1",  "IL8",  "CACNA1C",  "DSE",  "FAM84A",  "SLC16A2",  "COL6A2",  "FZD2",  "EFNB2",  "LOXL1-AS1",  
         "ARTN",  "STC1",  "LTBP1",  "MECOM",  "ZNF804A",  "RARRES1",  "RELN",  "TMEM200B",  "LAPTM5",  "MMP2",  "CAMK1D",  "PDGFC",  "PODXL",  "GCNT1",  "FGF2",  "CD96",  "CYR61",  "HIC1",  "MPP7",  "A4GALT",  "LINC00960",  "LAYN",  "WISP1",  "MX1",  "NPTX1",  "WNT5A-AS1",  "WNT5B",  "AC005789.11",  "FAM65C",  "IRF1",  "SAMD9L",  "KCNT2",  "FOXE1",  "NHS",  "LY6K",  "FOXC2",  "CRIM1",  "DOCK3",  "TNFRSF10A",  "DOCK2",  "CYP27C1",  "AMZ1",  "SECTM1",  "SOX9",  "FIBCD1",  "SH3RF2",  "DCN",  "SCARA3",  "BLK",  "TAGLN",  "PCOLCE",  "THSD4",  "KIAA1462",  "DKK1",  "FAM167A",  "PCDHGB4",  
         "SAMD12",  "B3GNT9",  "IFI44",  "DPYD",  "TLE4",  "IL20RB",  "FAM20C",  "PDGFB",  "SYT1",  "C12orf75",  "SCG5",  "EPHA2",  "ADAM8",  "AMPD3",  "LTBP2",  "MAGEC2",  "EPS8L2",  "RYR2",  "PCDHB3",  "NNMT",  "VASN",  "SPEG",  "LPHN2",  "PCDHGA11",  "PLAUR",  "ITGB4",  "ITGA2",  "CXCL3",  "ADAMTS4",  "PRDM1",  "STEAP2",  "BASP1",  "SPRED3",  "S100A16",  "TDRD9",  "MYLK",  "COL6A1",  "ANGPTL4",  "PCDH10",  "TFPI",  "PDZD2",  "RP11-575F12.2",  "SUSD1",  "HTR7",  "TPM2",  "CPE",  "FLNB",  "ZNF415",  "GPX1P1",  "SYBU",  "C1R",  "ERRFI1",  "IGFBP3",  "DPF3",  "ZNF185",  "CXCL1",  "PIP5KL1",  "GPR176",  "TRABD2A",  "NR2F1-AS1",  "IFI6",  "VCAN",  "SHC3",  "SCN2A",  "NUAK2",  "ITGA5",  "SPOCK1",  "IKZF2",  "DDAH1",  "BCL3",  "TFPI2",  "INHBE",  "RRAS",  "OSMR",  "LAMB3",  "EFEMP2",  "NMI",  "IFIT1")
geneSets <- GeneSet(genes, setName="Verfaillie_INV")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_malignant_subset_BT@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Verfaillie_INV<-getAUC(cells_AUC)
Verfaillie_INV<-t(Verfaillie_INV)
GC_malignant_subset_BT@meta.data<-cbind(GC_malignant_subset_BT@meta.data, Verfaillie_INV)
pdf("/Verfaillie_INV.pdf", width = 7.08, height = 5.8)
FeaturePlot(GC_malignant_subset_BT, features = "Verfaillie_INV",label = T)
dev.off()
pdf("/Verfaillie_INV_vln.pdf", width = 5, height = 4)
VlnPlot(GC_malignant_subset_BT, features = "Verfaillie_INV", pt.size = 0, group.by = "Malignant_clusters") + theme(plot.margin=unit(c(1,1,1.5,1.5),"cm"), legend.position = 'none')
dev.off()

##########################        Verfaillie_IMMUNE        ##########################
genes<-c("CLEC4C",  "GPAT2",  "LOC100133893",  "PTPRC",  "SIAH3",  "ABCB1",  "ABCB11",  "ABCC3",  "ABCD2",  "ABI3",  "ABI3BP",  "ACAP1",  "ACOXL",  "ACRBP",  "ACSL5",  "ACSL6",  "ACSM5",  "ACY3",  "ADAM28",  "ADAM6",  "ADAM8",  "ADAMDEC1",  "ADAP2",  "ADCY7",  "ADH1B",  "ADORA2A",  "ADORA3",  "ADRA1A",  "ADRA2A",  "AGAP2",  "AHSP",  "AICDA",  "AIF1",  "AIM2",  "AIPL1",  "AIRE",  "AKAP14",  "AKAP5",  "AKNA",  "ALOX5",  "ALOX5AP",  "ALPK2",  "AMICA1",  "AMPD1",  "AMPD3",  "ANK1",  "ANKRD22",  "ANKRD26P1",  "ANKRD29",  "ANKRD34B",  "ANKRD55",  "ANKRD58",  "ANO9",  "AOAH",  "AOX1",  "APBB1IP",  "APOB48R",  "APOBEC3A",  "APOBEC3D",  "APOBEC3G",  "APOBEC3H",  "APOF",  "APOL1",  "APOL3",  "APOL4",  "APOL6",  "ARHGAP15",  "ARHGAP25",  "ARHGAP27",  "ARHGAP30",  "ARHGAP9",  "ARHGDIB",  "ARL11",  "ARRDC5",  "ART4",  "ASCL2",  "ASGR2",  "ASMT",  "ATHL1",  "ATP1A4",  "ATP2A3",  "ATP8B4",  "B2M",  "BANK1",  "BASP1",  "BATF",  "BATF2",  "BCAS4",  "BCL11A",  "BCL11B",  "BCL2L14",  "BCORL2",  "BEND4",  "BET3L",  "BFSP2",  "BHLHA15",  "BIN2",  "BIRC3",  "BLK",  "BLNK",  "BMP10",  "BTG2",  "BTK",  "BTLA",  "BTN1A1",  "BTN3A1",  "BTN3A3",  "BTNL3",  "BTNL8",  "C10ORF105",  "C10ORF128",  "C10ORF50",  "C10ORF54",  "C11ORF21",  "C11ORF75",  "C12ORF42",  "C12ORF59",  "C12ORF63",  "C12ORF74",  "C12ORF77",  "C13ORF30",  "C14ORF115",  "C14ORF64",  "C14ORF73",  "C15ORF50",  "C15ORF53",  "C16ORF54",  "C17ORF28",  "C17ORF46",  "C17ORF60",  "C17ORF66",  "C17ORF87",  "C17ORF99",  "C19ORF38",  "C1ORF110",  "C1ORF162",  "C1ORF186",  "C1ORF200",  "C1ORF228",  "C1ORF54",  "C1QA",  "C1QB",  "C1QC",  "C1R",  "C1S",  "C2",  "C20ORF118",  "C20ORF141",  "C20ORF197",  "C21ORF96",  "C2ORF53",  "C2ORF85",  "C3",  "C3AR1",  "C3ORF51",  "C4A",  "C4BPB",  "C4ORF22",  "C4ORF50",  "C4ORF7",  "C5ORF20",  "C5ORF39",  "C5ORF40",  "C5ORF56",  "C6ORF97",  "C7",  "C7ORF72",  "C8ORF80",  "C9ORF128",  "C9ORF139",  "C9ORF144",  "CACNA1I",  "CACNG3",  "CAMK1D",  "CARD11",  "CARD17",  "CARD9",  "CASP10",  "CASP5",  "CASS4",  "CBFA2T3",  "CCDC109B",  "CCDC141",  "CCDC42",  "CCDC64",  "CCDC69",  "CCDC88B",  "CCDC88C",  "CCL1",  "CCL14",  "CCL16",  "CCL19",  "CCL2",  "CCL21",  "CCL22",  "CCL24",  "CCL25",  "CCL3",  "CCL3L1",  "CCL3L3",  "CCL4",  "CCL4L2",  "CCL5",  "CCL7",  "CCL8",  "CCND2",  "CCNI2",  "CCR1",  "CCR2",  "CCR3",  "CCR4",  "CCR5",  "CCR6",  "CCR7",  "CCR8",  "CCR9",  "CCRL2",  "CD14",  "CD160",  "CD163",  "CD180",  "CD19",  "CD1B",  "CD1C",  "CD1D",  "CD1E",  "CD2",  "CD200R1",  "CD200R1L",  "CD209",  "CD22",  "CD226",  "CD244",  "CD247",  "CD27",  "CD274",  "CD28",  "CD300A",  "CD300C",  "CD300E",  "CD300LB",  "CD300LF",  "CD33",  "CD37",  "CD38",  "CD3D",  "CD3E",  "CD3G",  "CD4",  "CD40",  "CD40LG",  "CD48",  "CD5",  "CD52",  "CD53",  "CD5L",  "CD6",  "CD69",  "CD7",  "CD70",  "CD72",  "CD74",  "CD79A",  "CD79B",  "CD80",  "CD84",  "CD86",  "CD8A",  "CD8B",  "CD96",  "CDC42SE2",  "CDH17",  "CDX1",  "CEACAM21",  "CEACAM22P",  "CEACAM3",  "CEACAM4",  "CEBPA",  "CEBPE",  "CECR1",  "CELA1",  "CFB",  "CFHR1",  "CFP",  "CHIT1",  "CHST15",  "CHST4",  "CIB3",  "CIITA",  "CLC",  "CLEC10A",  "CLEC12A",  "CLEC17A",  "CLEC4A",  "CLEC4D",  "CLEC4E",  "CLEC4M",  "CLEC6A",  "CLEC7A",  "CLEC9A",  "CLECL1",  "CLIC2",  "CLIC5",  "CLLU1OS",  "CLNK",  "CLU",  "CMAH",  "CMKLR1",  "CMPK2",  "CNR2",  "COL4A3",  "COL4A4",  "COLQ",  "CORO1A",  "CP",  "CPA5",  "CPA6",  "CPNE9",  "CPXCR1",  "CR1",  "CR1L",  "CR2",  "CREB3L3",  "CRIP3",  "CRTAM",  "CRYBB1",  "CSF1",  "CSF1R",  "CSF2RB",  "CSF3R",  "CSN3",  "CST2",  "CST7",  "CTSS",  "CTSW",  "CUX2",  "CX3CR1",  "CXCL10",  "CXCL11",  "CXCL12",  "CXCL13",  "CXCL16",  "CXCL9",  "CXCR2P1",  "CXCR3",  "CXCR4",  "CXCR5",  "CXCR6",  "CXORF21",  "CXORF50B",  "CYBB",  "CYLD",  "CYP2C8",  "CYSLTR1",  "CYTH4",  "CYTIP",  "DAO",  "DAPP1",  "DARC",  "DAZL",  "DBH",  "DCAF8L1",  "DCC",  "DDX4",  "DDX60",  "DEF6",  "DEFA4",  "DENND1C",  "DENND2D",  "DENND3",  "DERL3",  "DGKA",  "DHRS9",  "DIO1",  "DIO3",  "DIO3OS",  
         "DMRT3", "DMRTC1B",  "DNAH8",  "DNAJC5B",  "DNASE1L3",  "DNTT",  "DOCK11",  "DOCK2",  "DOCK8",  "DOK2",  "DOK3",  "DPEP1",  "DPEP2",  "DPPA4",  "DPT",  "DRAM1",  "DRD5",  "DTHD1",  "DTX1",  "DUPD1",  "DUSP2",  "DUSP26",  "EAF2",  "EBI3",  "ECEL1",  "EDAR",  "EGOT",  "EMB",  "EMR1",  "EMR4P",  "ENPP7",  "EOMES",  "EPHA1",  "EPSTI1",  "ERAP2",  "ERMN",  "ETV7",  "EVI2A",  "EVI2B",  "EVL",  "FAAH2",  "FAIM3",  "FAM107B",  "FAM113B",  "FAM129C",  "FAM151A",  "FAM153A",  "FAM153B",  "FAM157B",  "FAM159A",  "FAM163B",  "FAM166B",  "FAM170B",  "FAM177B",  "FAM179A",  "FAM18A",  "FAM26F",  "FAM46C",  "FAM49A",  "FAM55D",  "FAM65C",  "FAM92B",  "FASLG",  "FBP1",  "FCAMR",  "FCAR",  "FCER1G",  "FCER2",  "FCGR1A",  "FCGR1B",  "FCGR1C",  "FCGR3A",  "FCHO1",  "FCN1",  "FCRL1",  "FCRL2",  "FCRL3",  "FCRL4",  "FCRL5",  "FCRL6",  "FER1L4",  "FERMT3",  "FGD2",  "FGD3",  "FGL2",  "FGR",  "FLI1",  "FLJ40330",  "FLJ41941",  "FLJ43390",  "FLJ45983",  "FLT3",  "FLT3LG",  "FMNL1",  "FMO2",  "FMO6P",  "FOLR4",  "FOXD4L3",  "FOXP3",  "FPR1",  "FPR2",  "FPR3",  "FUCA1",  "FUT7",  "FXYD2",  "FYB",  "GAB3",  "GADD45G",  "GAGE12D",  "GAPT",  "GATA3",  "GBA3",  "GBP1",  "GBP2",  "GBP3",  "GBP4",  "GBP5",  "GBP7",  "GCET2",  "GCH1",  "GCNT1",  "GCNT4",  "GFI1",  "GGTA1",  "GHRL",  "GIMAP1",  "GIMAP2",  "GIMAP4",  "GIMAP5",  "GIMAP6",  "GIMAP7",  "GIMAP8",  "GJA10",  "GJD3",  "GLP2R",  "GLRX",  "GMFG",  "GNG8",  "GNGT2",  "GNLY",  "GOLGA2B",  "GP1BA",  "GP5",  "GPA33",  "GPBAR1",  "GPR114",  "GPR120",  "GPR132",  "GPR141",  "GPR142",  "GPR15",  "GPR160",  "GPR171",  "GPR174",  "GPR18",  "GPR182",  "GPR183",  "GPR25",  "GPR31",  "GPR35",  "GPR65",  "GPR68",  "GPR82",  "GPR84",  "GPRC5D",  "GRAP",  "GRAP2",  "GRAPL",  "GRIN3A",  "GSDMB",  "GTSF1",  "GTSF1L",  "GUCY2C",  "GVIN1",  "GZMA",  "GZMB",  "GZMH",  "GZMK",  "GZMM",  "H2-AA",  "H2-AB1",  "H2-BI",  "H2-D1",  "H2-DMA",  "H2-DMB1",  "H2-DMB2",  "H2-EB1",  "H2-EB2",  "H2-M1",  "H2-M10.1",  "H2-M10.2",  "H2-M10.3",  "H2-M10.4",  "H2-M10.5",  "H2-M10.6",  "H2-M11",  "H2-M2",  "H2-M3",  "H2-M5",  "H2-M9",  "H2-OA",  "H2-OB",  "H2-Q1",  "H2-Q10",  "H2-Q2",  "H2-Q4",  "H2-Q6",  "H2-Q7",  "H2-T10",  "H2-T11",  "H2-T22",  "H2-T23",  "H2-T24",  "H2-T3",  "H2-T9",  "HAMP",  "HAPLN3",  "HAVCR1",  "HAVCR2",  "HBD",  "HCG26",  "HCG27",  "HCK",  "HCLS1",  "HCP5",  "HCST",  "HEMGN",  "HIST1H1D",  "HIST1H2BI",  "HK3",  "HMGCLL1",  "HMHA1",  "HMHB1",  "HMSD",  "HP",  "HPD",  "HPR",  "HRASLS2",  "HSD11B1",  "HSF5",  "HSH2D",  "HTR3A",  "HTRA4",  "HVCN1",  "ICAM2",  "ICAM3",  "ICOS",  "IDO1",  "IDO2",  "IFI27",  "IFI30",  "IFI44L",  "IFIH1",  "IFIT3",  "IFNA21",  "IFNB1",  "IFNG",  "IGJ",  "IGLL1",  "IGSF6",  "IKZF1",  "IKZF2",  "IKZF3",  "IL10",  "IL10RA",  "IL12B",  "IL12RB1",  "IL15",  "IL15RA",  "IL17C",  "IL17REL",  "IL18",  "IL18BP",  "IL18R1",  "IL18RAP",  "IL2",  "IL21",  "IL21R",  "IL22RA2",  "IL23A",  "IL26",  "IL27",  "IL29",  "IL2RA",  "IL2RB",  "IL2RG",  "IL32",  "IL33",  "IL34",  "IL3RA", 
         "IL4",  "IL4I1",  "IL5RA",  "IL7",  "IL7R",  "IL9R",  "ILDR1",  "INMT",  "INPP5D",  "INSL3",  "INSM1",  "IPCEF1",  "IQCF1",  "IRF1",  "IRF5",  "IRF8",  "IRGM",  "ISG15",  "ITGAD",  "ITGAL",  "ITGAM",  "ITGAX",  "ITGB2",  "ITGB7",  "ITIH1",  "ITIH3",  "ITIH4",  "ITK",  "ITM2A",  "JAK2",  "JAK3",  "JAKMIP1",  "JSRP1",  "KBTBD8",  "KCNA3",  "KCNJ1",  "KCNJ5",  "KCNK12",  "KCNK13",  "KCNN3",  "KEL",  "KIAA0125",  "KIAA0748",  "KIAA1324",  "KIAA2022",  "KIF19",  "KIF21B",  "KIR2DL1",  "KIR2DL3",  "KIR2DL4",  "KIR2DS4",  "KIR3DL1",  "KIR3DL2",  "KIR3DL3",  "KIR3DX1",  "KLHDC7B",  "KLHL33",  "KLHL6",  "KLRB1",  "KLRC1",  "KLRC2",  "KLRC3",  "KLRC4",  "KLRD1",  "KLRF1",  "KLRG1",  "KLRK1",  "KMO",  "KRT72",  "KSR2",  "LAG3",  "LAIR1",  "LAIR2",  "LAMP3",  "LAP3",  "LAPTM5",  "LAT",  "LAT2",  "LAX1",  "LBP",  "LCK",  "LCN8",  "LCP1",  "LCP2",  "LGALS14",  "LGALS2",  "LGALS9",  "LGALS9B",  "LGALS9C",  "LILRA1",  "LILRA2",  "LILRA3",  "LILRA4",  "LILRA5",  "LILRA6",  "LILRB1",  "LILRB2",  "LILRB3",  "LILRB4",  "LILRB5",  "LILRP2",  "LIM2",  "LIMD2",  "LIME1",  "LINGO3",  "LINGO4",  "LMO2",  "LOC100124692",  "LOC100128164",  "LOC100128542",  "LOC100129066",  "LOC100130933",  "LOC100188949",  "LOC100192379",  "LOC100233209",  "LOC100272216",  "LOC145820",  "LOC254559",  "LOC283070",  "LOC283314",  "LOC283663",  "LOC283999",  "LOC284551",  "LOC284749",  "LOC285740",  "LOC285780",  "LOC285796",  "LOC400696",  "LOC400759",  "LOC440461",  "LOC440896",  "LOC606724",  "LOC647121",  "LOC729467",  "LOC96610",  "LPAR5",  "LPXN",  "LRMP",  "LRP2",  "LRRC25",  "LRRC55",  "LSP1",  "LST1",  "LTA",  "LTB",  "LTC4S",  "LY75",  "LY86",  "LY9",  "LYL1",  "LYN",  "LYZ",  "MACC1",  "MADCAM1",  "MAN1A1",  "MAP3K8",  "MAP4K1",  "MARCO",  "MATK",  "MDS2",  "MEF2B",  "MEFV",  "MEI1",  "MEP1A",  "MFNG",  "MGC29506",  "MIAT",  "MIR155HG",  "MIXL1",  "MLKL",  "MMP12",  "MMP25",  "MMRN1",  "MNDA",  "MOGAT2",  "MPEG1",  "MRC1",  "MS4A1",  "MS4A4A",  "MS4A6A",  "MSR1",  "MSX2P1",  "MTTP",  "MUSK",  "MYB",  "MYBPC2",  "MYBPC3",  "MYO1A",  "MYO1F",  "MYO1G",  "MYO7A",  "NAP1L2",  "NAPSA",  "NAPSB",  "NAT8B",  "NCF1",  "NCF1B",  "NCF1C",  "NCF2",  "NCF4",  "NCKAP1L",  "NCR1",  "NCR2",  "NCR3",  "NCRNA00158",  "NCRNA00161",  "NETO1",  "NEURL3",  "NEUROG3",  "NFAM1",  "NKG7",  "NKX6-3",  "NLRC3",  "NLRC4",  "NLRC5",  "NLRP1",  "NLRP3",  "NLRP6",  "NLRP7",  "NMI",  "NOD2",  "NOL4",  "NPHS1",  "NPY1R",  "NPY5R",  "NTNG2",  "NUAK2",  "OAS1",  "OAS2",  "OASL",  "OCM",  "ODF3B",  "OLR1",  "OR10G2",  "OR13A1",  "OR3A3",  "OR4C6",  "OR52N2",  "OR52N4",  "OR56B1",  "OR5C1",  "OSCAR",  "OTOF",  "P2RX1",  "P2RX2",  "P2RX5",  "P2RY10",  "P2RY12",  "P2RY13",  "P2RY14",  "P2RY6",  "P2RY8",  "PACSIN1",  "PADI4",  "PARP12",  "PARP14",  "PARP15",  "PARVG",  "PATL2",  "PAX5",  "PCDH11Y",  "PCDHAC1",  "PDCD1",  "PDCD1LG2",  "PDE6H",  "PECAM1",  "PIK3AP1",  "PIK3CG",  "PIK3IP1",  "PIK3R5",  "PIK3R6",  "PILRA",  "PIM1",  "PIM2",  "PKD2L1",  "PKHD1L1",  "PLA2G2D",  "PLA2G7",  
         "PLAC8", "PLCB2",  "PLCG2",  "PLCH2",  "PLCL2",  "PLCXD2",  "PLD4",  "PLD5",  "PLEK",  "PLEKHG7",  "PM20D1",  "PMCH",  "PNOC",  "POM121L9P",  "POU2AF1",  "POU2F2",  "PPP1R16B",  "PPP1R2P9",  "PRAM1",  "PRDM1",  "PRDM8",  "PRF1",  "PRKCB",  "PRKCH",  "PRKCQ",  "PROC",  "PRODH2",  "PROKR2",  "PRR5L",  "PRSS1",  "PSD4",  "PSG2",  "PSG5",  "PSMB10",  "PSMB8",  "PSMB9",  "PSME2",  "PSTPIP1",  "PSTPIP2",  "PTAFR",  "PTGDR",  "PTGER2",  "PTGER4",  "PTK2B",  "PTPN22",  "PTPN6",  "PTPN7",  "PTPRCAP",  "PTPRN2",  "PTPRO",  "PTPRQ",  "PVALB",  "PVRIG",  "PYHIN1",  "RAB19",  "RAB37",  "RAC2",  "RARRES1",  "RARRES3",  "RASAL3",  "RASGEF1A",  "RASGEF1B",  "RASGRP1",  "RASGRP2",  "RASSF4",  "RASSF5",  "RASSF6",  "RBP5",  "RCSD1",  "REC8",  "REG4",  "RFTN1",  "RGL4",  "RGS13",  "RGS18",  "RGS7BP",  "RGS9",  "RHOF",  "RHOH",  "RIC3",  "RIPK3",  "RLTPR",  "RNASE2",  "RNASE6",  "RNASET2",  "RSAD2",  "RSPO1",  "RSPO3",  "RUFY4",  "RUNDC2C",  "S100Z",  "S1PR1",  "S1PR4",  "SAMD3",  "SAMD7",  "SAMD9L",  "SAMHD1",  "SAMSN1",  "SASH3",  "SCN3A",  "SCT",  "SDS",  "SECTM1",  "SEL1L3",  "SELL",  "SELP",  "SELPLG",  "SEMA4A",  "SEMA4D",  "SERPINA1",  "SERPINA9",  "SERPINF2",  "SERPING1",  "SH2D1A",  "SH2D1B",  "SH2D2A",  "SH2D3C",  "SHD",  "SHISA3",  "SIDT1",  "SIGLEC1",  "SIGLEC10",  "SIGLEC11",  "SIGLEC14",  "SIGLEC16",  "SIGLEC5",  "SIGLEC6",  "SIGLEC7",  "SIGLEC8",  "SIGLEC9",  "SIRPB2",  "SIRPD",  "SIRPG",  "SIT1",  "SKAP1",  "SLA",  "SLA2",  "SLAMF1",  "SLAMF6",  "SLAMF7",  "SLAMF8",  "SLC12A3",  "SLC15A3",  "SLC1A2",  "SLC22A2",  "SLC22A3",  "SLC27A2",  "SLC29A3",  "SLC2A5",  "SLC32A1",  "SLC40A1",  "SLC4A1",  "SLC5A5",  "SLC6A7",  "SLC7A7",  "SLC8A1",  "SLC9A9",  "SLCO2B1",  "SLFN12L",  "SLFN13",  "SLFN14",  "SMTNL1",  "SNAI3",  "SNAP91",  "SNX20",  "SOCS1",  "SOX30",  "SP140",  "SPATC1",  "SPI1",  "SPIB",  "SPIC",  "SPINK2",  "SPN",  "SPNS3",  "SPOCK2",  "SPTA1",  "SQRDL",  "SRGN",  "SSTR2",  "SSTR3",  "SSX2",  "ST8SIA4",  "STAB2",  "STAC3",  "STAP1",  "STAT1",  "STAT4",  "STK17B",  "STX11",  "STYK1",  "SUCNR1",  "SULT1B1",  "SUMO1P1",  "SUSD3",  "SYCE1",  "SYK",  "SYT15",  "SYTL1",  "SYTL3",  "T",  "TAGAP",  "TAP1",  "TARP",  "TAS1R3",  "TBC1D10C",  "TBX21",  "TCF7",  "TCL1A",  "TCL1B",  "TDGF1",  "TDGF3",  "TDRD1",  "TECRL",  "TFEC",  "THEMIS",  "TIFAB",  "TIGIT",  "TIMD4",  "TLL1",  "TLR10",  "TLR7",  "TLR8",  "TLR9",  "TMC8",  "TMEM130",  "TMEM132C",  "TMEM149",  "TMEM150B",  "TMEM155",  "TMEM156",  "TMEM171",  "TMEM176A",  "TMEM176B",  "TMEM37",  "TMEM90B",  "TMIGD2",  "TMPRSS3",  "TMSL3",  "TNF",  "TNFAIP2",  "TNFAIP8",  "TNFAIP8L2",  "TNFRSF13B",  "TNFRSF13C",  "TNFRSF17",  "TNFRSF18",  "TNFRSF1B",  "TNFRSF25",  "TNFRSF9",  "TNFSF10",  "TNFSF11",  "TNFSF12-TNFSF13",  "TNFSF13B",  "TNFSF14",  "TNFSF8",  "TNIP3",  "TNNI2",  "TOX",  "TOX2",  "TRAF3IP3",  "TRANK1",  "TRAT1",  "TREM2",  "TREML1",  "TREML2",  "TREML4",  "TRIM22",  "TRIM69",  "TRPC2",  "TSHR",  "TSPAN32",  "TTBK1",  "TTC16",  "TTC24",  "TTC7A",  "TTC9",  "TXK",  "TXNDC3",  "TYMP",  "TYROBP",  "UBA7",  "UBASH3A",  "UBD",  "UBE2L6",  "UCP2",  "UGT2B15",  "UNC13D",  "UNC45B",  "UNC5A",  "UPK3A",  "UTS2",  "UTS2R",  "VAMP5",  "VAV1",  "VCAM1",  "VIPR1",  "VIPR2",  "VNN1",  "VNN2",  "VOPP1",  "VPREB1",  "VPREB3",  "VSIG1",  "VSIG4",  "VSTM1",  "VSX2",  "VWA5B1",  "WARS",  "WAS",  "WDFY4",  "WDR49",  "WNT1",  "WNT10A",  "XAF1",  "XCL1",  "XCL2",  "XCR1",  "XIRP1",  "XKR4",  "ZAP70",  "ZBED2",  "ZBP1",  "ZBTB32",  "ZC3H12D",  "ZMYND15",  "ZNF385D",  "ZNF683",  "ZNF80",  "ZNF804A",  "ZNF831",  "ZPBP2")
geneSets <- GeneSet(genes, setName="Verfaillie_IMMUNE")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_malignant_subset_BT@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Verfaillie_IMMUNE<-getAUC(cells_AUC)
Verfaillie_IMMUNE<-t(Verfaillie_IMMUNE)
GC_malignant_subset_BT@meta.data<-cbind(GC_malignant_subset_BT@meta.data, Verfaillie_IMMUNE)
pdf("/Verfaillie_IMMUNE.pdf", width = 7.08, height = 5.8)
FeaturePlot(GC_malignant_subset_BT, features = "Verfaillie_IMMUNE",label = T)
dev.off()
pdf("/Verfaillie_IMMUNE_vln.pdf", width = 5, height = 4)
VlnPlot(GC_malignant_subset_BT, features = "Verfaillie_IMMUNE", pt.size = 0, group.by = "Malignant_clusters") + theme(plot.margin=unit(c(1,1,1.5,1.5),"cm"), legend.position = 'none')
dev.off()

##########################        Wouters_Intermediate        ##########################
gmtFile <- paste("/Wouters_NTERMEDIATE.gmt")
geneSets <- getGmt(gmtFile)
geneSets
#cells_rankings <- AUCell_buildRankings(GC_malignant_subset_BT@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Wouters_Intermediate<-getAUC(cells_AUC)
Wouters_Intermediate<-t(Wouters_Intermediate)
GC_malignant_subset_BT@meta.data<-cbind(GC_malignant_subset_BT@meta.data, Wouters_Intermediate)
pdf("/Wouters_Intermediate.pdf", width = 7.08, height = 5.8)
FeaturePlot(GC_malignant_subset_BT, features = "Wouters_Intermediate",label = T)
dev.off()
pdf("/Wouters_Intermediate_vln.pdf", width = 5, height = 4)
VlnPlot(GC_malignant_subset_BT, features = "Wouters_Intermediate", pt.size = 0, group.by = "Malignant_clusters") + theme(plot.margin=unit(c(1,1,1.5,1.5),"cm"), legend.position = 'none')
dev.off()

##########################        Wouters_Melanocytic        ##########################
gmtFile <- paste("/Wouters_melanocytic.gmt")
geneSets <- getGmt(gmtFile)
geneSets
#cells_rankings <- AUCell_buildRankings(GC_malignant_subset_BT@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Wouters_Melanocytic<-getAUC(cells_AUC)
Wouters_Melanocytic<-t(Wouters_Melanocytic)
GC_malignant_subset_BT@meta.data<-cbind(GC_malignant_subset_BT@meta.data, Wouters_Melanocytic)
pdf("/Wouters_Melanocytic.pdf", width = 7.08, height = 5.8)
FeaturePlot(GC_malignant_subset_BT, features = "Wouters_Melanocytic",label = T)
dev.off()
pdf("/Wouters_Melanocytic_vln.pdf", width = 5, height = 4)
VlnPlot(GC_malignant_subset_BT, features = "Wouters_Melanocytic", pt.size = 0, group.by = "Malignant_clusters") + theme(plot.margin=unit(c(1,1,1.5,1.5),"cm"), legend.position = 'none')
dev.off()

##########################        Wouters_MESlike        ##########################
gmtFile <- paste("/Wouters_MESlike.gmt")
geneSets <- getGmt(gmtFile)
geneSets
#cells_rankings <- AUCell_buildRankings(GC_malignant_subset_BT@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Wouters_MESlike<-getAUC(cells_AUC)
Wouters_MESlike<-t(Wouters_MESlike)
GC_malignant_subset_BT@meta.data<-cbind(GC_malignant_subset_BT@meta.data, Wouters_MESlike)
pdf("/Wouters_MESlikec.pdf", width = 7.08, height = 5.8)
FeaturePlot(GC_malignant_subset_BT, features = "Wouters_MESlike",label = T)
dev.off()
pdf("/Wouters_MESlike_vln.pdf", width = 5, height = 4)
VlnPlot(GC_malignant_subset_BT, features = "Wouters_MESlike", pt.size = 0, group.by = "Malignant_clusters") + theme(plot.margin=unit(c(1,1,1.5,1.5),"cm"), legend.position = 'none')
dev.off()


####################################heatmap#########################################################
reglist <- c("Malignant_clusters","Hoek_INV","Hoek_PRO","Tirosh_AXL","Tirosh_MITF","Tsoi_Melanocytic","Tsoi_Transit","Tsoi_undiff", "Tsoi_NClike","Verfaillie_IMMUNE","Verfaillie_INV","Verfaillie_PRO","Wouters_Intermediate","Wouters_Melanocytic","Wouters_MESlike","BARON_Melanocytic","BARON_NClike","BARON_STRESS","Rambow_Hyperdiff","Rambow_MESlike","Rambow_NCSC","Rambow_SMC")
matr <- FetchData(GC_malignant_subset_BT, vars = reglist, cells = NULL) # %>%
#apply(2, function(x) (x - min(x)) / (max(x) - min(x))) # %>%
#t
library(scales)
matr <- matr[ , !duplicated(colnames(matr))]
matr <- matr %>% group_by(Malignant_clusters) %>% summarise_all(mean)
#
Malignant_clusters <- levels(GC_malignant_subset_BT$Malignant_clusters)
Malignant_clusters
color_1 <- setNames(hue_pal()(11), Malignant_clusters)
color_list <- list(color_1 = Malignant_clusters) 
color_list
#names(color_list) <- list(Malignant_clusters_old)
Malignant_clusters <- matr$Malignant_clusters
ha = ComplexHeatmap::HeatmapAnnotation(Malignant_clusters = Malignant_clusters,
                       annotation_name_side = "left",
                       annotation_name_rot = 180,
                       show_annotation_name = FALSE
                       #col = color_list
)

#head(matr)
#matr <-apply(matr, 1, as.numeric)
head(matr)

matr$Malignant_clusters <- NULL
mat_num <- matrix(as.numeric(as.character((matr),    # Convert to numeric matrix
                                          ncol = ncol(matr))))
matr <- t(matr)
mat_scaled = t(apply(matr, 1, scale))
ht<-ComplexHeatmap::Heatmap(mat_scaled, 
            cluster_rows = TRUE, 
            cluster_columns = FALSE,
            show_column_names = FALSE,
            top_annotation = ha,
            #col = color_list,
            column_split = Malignant_clusters)
pdf("/scores_heatmap_av.pdf", width = 10, height = 6)
draw(ht)
dev.off()



#################  Human (new) signatures

##########################        Neural_crest_like        ##########################
genes<-c("DCSTAMP",  "COL8A1",  "NGFR",  "AXL",  "MRC2",  "BMP8B",  "ECE1",  "NTM",  "TNC",  "PMEPA1",  "FREM2",  "ITGB8",  "COL4A1",  "SLC20A1",  "LOXL3",  "OLFML3",  "COL5A2",  "ITGA3",  "IL1RAP",  "ITIH5",  "IGFBP3",  "COL12A1",  "THBS2",  "KCNN4",  "COL4A2",  "LYPD1",  "CALU",  "KCNMA1",  "FLRT3",  "ITGA6",  "ITGB3",  "ACTN4",  "MAP1B",  "TNFRSF19",  "PTN",  "P4HB",  "SCN9A",  "ADAMTS1",  "FNDC3B",  "ITGB1",  "LTBP1",  "NT5E",  "SERPINE2",  "PCDH9",  "COL9A3",  "DCBLD2",  "CHST6",  "CADM1",  "HSP90B1",  "SEMA3B",  "TGFBI",  "DAG1",  "TIMP1",  "CDCP1",  "TIMP3",  "PRNP",  "MICALL2",  "EMILIN1",  "PDIA4",  "PPFIBP1",  "COL15A1",  "CREB5",  "IGFBP5",  "FN1",  "PDIA6",  "ITGA1",  "CDH2",  "SPARC",  "RAB31",  "PDIA3",  "PLEC",  "COL18A1",  "SEMA3A",  "AMOTL1",  "COL1A2",  "HYOU1",  "PLAT",  "TANC2",  "PDGFA",  "SDC3",  "DPY19L1",  "FAM3C",  "CD109",  "ALDH1A3",  "JAG1",  "PLAUR",  "RUNX1",  "COL6A2",  "SKI",  "COL6A1",  "SEC24D",  "TNFRSF12A",  "SH3PXD2A",  "SCARB2",  "APLP2",  "RTN4",  "PLOD2",  "ECM1",  "CREB3L2")
geneSets <- GeneSet(genes, setName="Neural_crest_like")
geneSets
cells_rankings <- AUCell_buildRankings(GC_malignant_subset_BT@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Neural_crest_like<-getAUC(cells_AUC)
Neural_crest_like<-t(Neural_crest_like)
GC_malignant_subset_BT@meta.data<-cbind(GC_malignant_subset_BT@meta.data, Neural_crest_like)
##########################        Interferon_responsive        ##########################
genes<-c("DPP6",  "PKHD1",  "LINC00323",  "PCSK2",  "IFIT1",  "OAS2",  "PALMD",  "TFAP2B",  "CTXND1",  "APCDD1",  "LGI3",  "EIF1AY",  "ANGPT2",  "MX1",  "SLC1A4",  "LY6E",  "GNG11",  "ISG15",  "IFIT3",  "OAS1",  "RAMP1",  "ABHD2",  "HR",  "IFI27",  "SEMA3D",  "IFI44L",  "SFTPC",  "SCUBE2",  "XAF1",  "GMPR",  "CMPK2",  "KCNE4",  "HERC6",  "CD9",  "SOX9",  "MAL",  "ADCY2",  "IFI44",  "NFATC2",  "MX2",  "USP18",  "TGFBI",  "IFI6",  "IL16",  "BST2",  "SLC9A3",  "CAV1",  "PRELP",  "DDX3Y",  "FAM69C",  "SAMD9L",  "TYRP1",  "C10orf90",  "CST3",  "CAVIN1",  "WDR63",  "PML",  "C9orf3",  "QPCT",  "DSTYK",  "MT1G",  "IFIT2",  "MAPK8IP2",  "PERP",  "SOD1",  "LY6E-DT",  "OAS3",  "ATP6V0E2",  "CYB561A3",  "SP110",  "EMP3",  "SLC25A4",  "RPS4Y1",  "EZR",  "PLXNC1",  "CTSH",  "SLCO4A1",  "COPZ2",  "VEPH1",  "SORBS2",  "ACSL1",  "VAMP5",  "SYNM",  "CTSB",  "CAV2",  "CSPG4",  "SAP25",  "RSAD2",  "SPATS2L",  "CYCS",  "NEO1",  "LGALS3BP",  "ANXA2",  "HTRA1",  "IRF7",  "BMP2K",  "CISD1",  "GJB1",  "CSTB",  "ANKRD9")
geneSets <- GeneSet(genes, setName="Interferon_responsive")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_malignant_subset_BT@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Interferon_responsive<-getAUC(cells_AUC)
Interferon_responsive<-t(Interferon_responsive)
GC_malignant_subset_BT@meta.data<-cbind(GC_malignant_subset_BT@meta.data, Interferon_responsive)
##########################        Melanocytic        ##########################
genes<-c("RPS7",  "RPS18",  "RPS23",  "RPL7A",  "RPL12",  "RPS6",  "RPS13",  "RPS14",  "RPS3",  "RPL41",  "RPS15A",  "RPL9",  "RPS8",  "RPL8",  "NACA",  "RPL32",  "RPL30",  "RPL37",  "RPL6",  "RPL19",  "RPL35",  "RPL29",  "RPLP1",  "RPL18A",  "RPL5",  "RPS2",  "RPS3A",  "RPS15",  "RPL10A",  "RPS24",  "RPL35A",  "RACK1",  "RPL10",  "RPS27A",  "RPL18",  "RPL3",  "RPL14",  "RPS19",  "RPS4X",  "RPL24",  "RPL36",  "RPS9",  "RPLP2",  "RPL34",  "RPL15",  "RPS5",  "RPL11",  "RPL37A",  "RPS28",  "FAU",  "RPL28",  "RPS16",  "RPL39",  "RPL7",  "BTF3",  "RPSA",  "RPS27",  "RPLP0",  "RPL13",  "ATP5MC2",  "RPL31",  "RPL21",  "RPL23A",  "UBA52",  "ZFAS1",  "RPS25",  "RPS12",  "RPS21",  "EEF1A1",  "RPL36A",  "RPL4",  "RPL22",  "EIF3H",  "EIF3E",  "PFDN5",  "TPT1",  "EEF1B2",  "RPL27",  "RPL26",  "UQCRH",  "SLC25A3",  "EIF3F",  "RPL36AL",  "NPM1",  "RPS26",  "COX7C",  "FTH1",  "RPL17",  "UQCRB",  "EEF2",  "EIF3K",  "PABPC1",  "HNRNPA1",  "COX5A",  "SLC25A5",  "YBX1",  "TXN",  "FTL",  "RPL38",  "CLTA")
geneSets <- GeneSet(genes, setName="Melanocytic")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_malignant_subset_BT@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Melanocytic<-getAUC(cells_AUC)
Melanocytic<-t(Melanocytic)
GC_malignant_subset_BT@meta.data<-cbind(GC_malignant_subset_BT@meta.data, Melanocytic)
##########################        Mesenchymal        ##########################
genes<-c("PDGFRB",  "NRP1",  "COL5A1",  "COL1A1",  "COL3A1",  "LUM",  "THY1",  "TCF4",  "ID3",  "PMEPA1",  "MRC2",  "DCN",  "COL1A2",  "COL4A1",  "MXRA8",  "TMSB4X",  "COL6A2",  "COL6A3",  "COL5A2",  "COL4A2",  "COL6A1",  "AQP1",  "FBLN1",  "IGFBP7",  "COL12A1",  "BGN",  "PCOLCE",  "VCAN",  "ANTXR2",  "PBX1",  "HTRA1",  "C1R",  "SELENOM",  "FBN1",  "SPARC",  "MARCKS",  "FN1",  "COL15A1",  "PALLD",  "CALD1",  "ITGA1",  "POSTN",  "CFH",  "PTN",  "ITGB1",  "H3F3B",  "TPM2",  "TIMP1",  "TMEM47",  "PRRX1",  "PTMA",  "C1S",  "COL18A1",  "TUBA1A",  "TPM1",  "IGFBP2",  "AHNAK",  "CCDC80",  "TMSB10",  "MMP2",  "IGFBP4",  "CD9",  "PGRMC1",  "NR2F2",  "TPM4",  "CNN3",  "S100A4",  "LAMB1",  "EMP2",  "HMGB1",  "GSN",  "CTHRC1",  "TGFBR2",  "TNFRSF12A",  "FGFR1",  "S100A6",  "IFITM2",  "LGALS1",  "VCL",  "LAMC1",  "TGFBI",  "MGP",  "PTMS",  "ANXA2",  "ANXA1",  "LAPTM4A",  "GNG11",  "ITM2B",  "TIMP3",  "LHFPL6",  "H3F3A",  "FERMT2",  "DEK",  "FSTL1",  "ARGLU1",  "RAP1B",  "UTRN",  "MYL9",  "CAVIN1",  "NDRG2")
geneSets <- GeneSet(genes, setName="Mesenchymal")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_malignant_subset_BT@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Mesenchymal<-getAUC(cells_AUC)
Mesenchymal<-t(Mesenchymal)
GC_malignant_subset_BT@meta.data<-cbind(GC_malignant_subset_BT@meta.data, Mesenchymal)
##########################        Mitochondrial        ##########################
genes<-c("MGP",  "APOE",  "MT-ATP8",  "MT-ND6",  "DST",  "APP",  "PCDH7",  "NFIA",  "KCNJ10",  "RNASE1",  "MT-ND5",  "PTPRZ1",  "MOXD1",  "DDX17",  "NES",  "SPP1",  "MTRNR2L12",  "TBX2",  "XIST",  "ZEB2",  "MT-ND4L",  "ZBTB20",  "MT-ND3",  "RBMS3",  "AEBP1",  "Sep-07",  "CLTC",  "MPZL1",  "NBL1",  "TF",  "NLGN1",  "PLP1",  "MT-CO2",  "NORAD",  "MEGF10",  "ERBB3",  "TAOK1",  "MT-CO1",  "CEP170",  "FSTL1",  "DDX5",  "AKAP9",  "SNX22",  "IGSF3",  "HP1BP3",  "GPR37",  "LGI4",  "GLG1",  "ADGRG1",  "DCT",  "TRIM2",  "MTRNR2L8",  "ATRX",  "LPL",  "ITM2B",  "GPNMB",  "TEAD1",  "CTNNB1",  "ATP1A1",  "KIF13A",  "CHL1",  "CLIC4",  "CCDC50",  "SLC44A1",  "ALDOC",  "ENAH",  "PLAT",  "FRMD6",  "FADS2",  "PLCG2",  "APOC1",  "TSC22D1",  "SESN3",  "ZFHX3",  "ROBO1",  "BICD1",  "ITGAV",  "SYT11",  "ETV1",  "PIGM",  "FCRLA",  "PRKDC",  "PFKFB2",  "CDH19",  "TBL1XR1",  "SOX4",  "AC009041.2",  "AP1S2",  "SORT1",  "SPRY1",  "PTPN13",  "LMO4",  "UBL3",  "TNFRSF21",  "KCTD12",  "CC2D1A",  "GPRC5B",  "NMT1",  "BCAN",  "IGSF8")
geneSets <- GeneSet(genes, setName="Mitochondrial")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_malignant_subset_BT@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Mitochondrial<-getAUC(cells_AUC)
Mitochondrial<-t(Mitochondrial)
GC_malignant_subset_BT@meta.data<-cbind(GC_malignant_subset_BT@meta.data, Mitochondrial)
##########################        Antigen_presenting        ##########################
genes<-c("HLA-DRA",  "HLA-DPA1",  "HLA-DPB1",  "GBP1",  "HLA-DRB5",  "HLA-DQA1",  "HLA-B",  "HLA-DMA",  "HLA-DQB1",  "HLA-E",  "CIITA",  "HLA-DMB",  "TAP1",  "HLA-DRB1",  "HLA-C",  "RARRES3",  "PSMB9",  "STAT1",  "IL18BP",  "APOL6",  "B2M",  "GBP4",  "HLA-F",  "GBP2",  "APOL1",  "CD74",  "HLA-A",  "IRF1",  "WARS",  "IFIH1",  "TAPBP",  "LRP2",  "TAP2",  "PSMB8",  "PARP14",  "SERPING1",  "TRIM22",  "UBE2L6",  "ISG15",  "HAPLN3",  "PARP9",  "GBP3",  "XAF1",  "IFI35",  "NLRC5",  "VAMP5",  "CTSS",  "DTX3L",  "PSME1",  "SAMD9L",  "TNFSF13B",  "LAP3",  "PSME2",  "RSAD2",  "EPSTI1",  "ZNFX1",  "BTN3A1",  "PARP12",  "NMI",  "BIRC3",  "PDLIM4",  "OAS1",  "IFI6",  "C1R",  "CASP1",  "TNC",  "BST2",  "IFITM1",  "IFI44L",  "SAMHD1",  "BTN3A2",  "TMEM140",  "DDX60",  "BTN3A3",  "IRF7",  "GSDMD",  "OAS3",  "TRIM69",  "FRMD4A",  "A2M",  "IFITM3",  "TRIM56",  "RTP4",  "CTSO",  "STAT2",  "IFIT3",  "SP100",  "IFI16",  "NAMPT",  "OAS2",  "APOL2",  "CD47",  "CLEC2B",  "CAPG",  "HERC6",  "HELZ2",  "NUB1",  "MIA",  "GPR155",  "S100A10")
geneSets <- GeneSet(genes, setName="Antigen_presenting")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_malignant_subset_BT@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Antigen_presenting<-getAUC(cells_AUC)
Antigen_presenting<-t(Antigen_presenting)
GC_malignant_subset_BT@meta.data<-cbind(GC_malignant_subset_BT@meta.data, Antigen_presenting)
##########################        Mitotic        ##########################
genes<-c("E2F2",  "ESCO2",  "PKMYT1",  "ASF1B",  "RRM2",  "KIFC1",  "CLSPN",  "FAM111B",  "SPC25",  "EXO1",  "PBK",  "MYBL2",  "RAD51AP1",  "DTL",  "NCAPG",  "NCAPH",  "CDK1",  "CDCA5",  "NEIL3",  "TK1",  "CDC6",  "MCM10",  "POLQ",  "SKA3",  "NUF2",  "UHRF1",  "MKI67",  "ANLN",  "HJURP",  "MELK",  "NDC80",  "NUSAP1",  "TCF19",  "AURKB",  "CCNA2",  "FOXM1",  "CEP55",  "E2F1",  "KNL1",  "CDC25C",  "SPC24",  "SGO1",  "PLK4",  "CCNE2",  "FANCD2",  "DEPDC1",  "DIAPH3",  "CENPU",  "TOP2A",  "KIF15",  "FANCI",  "UBE2C",  "KIF2C",  "FBXO43",  "PCLAF",  "KIF18B",  "E2F7",  "WDR76",  "TACC3",  "MND1",  "CKAP2L",  "ASPM",  "BUB1B",  "KIF11",  "CENPM",  "CENPK",  "BLM",  "SHCBP1",  "RAD54L",  "IQGAP3",  "ZWINT",  "LMNB1",  "PRC1",  "SKA1",  "CDCA3",  "BRIP1",  "DLGAP5",  "BIRC5",  "ZNF367",  "OIP5",  "KIF14",  "MYBL1",  "CDCA2",  "CDCA8",  "CDC45",  "CHAF1A",  "BRCA2",  "TPX2",  "GTSE1",  "GINS2",  "HELLS",  "KIF4A",  "TRIP13",  "ATAD2",  "NCAPG2",  "BARD1",  "XRCC2",  "UBE2T",  "KNTC1",  "SPAG5")
geneSets <- GeneSet(genes, setName="Mitotic")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_malignant_subset_BT@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Mitotic<-getAUC(cells_AUC)
Mitotic<-t(Mitotic)
GC_malignant_subset_BT@meta.data<-cbind(GC_malignant_subset_BT@meta.data, Mitotic)
##########################        Stress_paraspeckles        ##########################
genes<-c("MALAT1",  "NEAT1",  "POLR2J3.1",  "GABPB1-AS1",  "SON",  "N4BP2L2",  "RBM39",  "PLEKHA5",  "FUS",  "HPS4",  "RSRP1",  "MDM4",  "RNMT",  "KCNQ1OT1",  "HNRNPU",  "SFPQ",  "CCNL1",  "PPP1R15A",  "ZNF704",  "VMP1",  "TCF25",  "ATP6V0A1",  "THUMPD3-AS1",  "TAF1D",  "SPG7",  "MIR3142HG",  "NKTR",  "ARID1B",  "SRSF10",  "ARGLU1",  "TRPM1",  "ATM",  "SMG1",  "CHD9",  "AC020916.1",  "ANKRD28",  "MICAL3",  "PNISR",  "PNN",  "IVNS1ABP",  "JUN",  "NHLRC3",  "GON4L",  "SF1",  "MCOLN3",  "RGS12",  "ERICH1",  "SLC25A37",  "INTS6",  "HMCN1",  "IGF2BP2",  "KANSL1",  "MYO5A",  "LUC7L3",  "MAML2",  "WSB1",  "SLC16A1-AS1",  "LYST",  "PAX3",  "HSPA1A",  "CALD1",  "RBM25",  "FMN1",  "IFI44L",  "CCNL2")
geneSets <- GeneSet(genes, setName="Stress_paraspeckles")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_malignant_subset_BT@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Stress_paraspeckles<-getAUC(cells_AUC)
Stress_paraspeckles<-t(Stress_paraspeckles)
GC_malignant_subset_BT@meta.data<-cbind(GC_malignant_subset_BT@meta.data, Stress_paraspeckles)
##########################        Stress_hypoxia        ##########################
genes<-c("VEGFA",  "NDUFA4L2",  "FLNB",  "NXPH4",  "SLC16A3",  "ENO2",  "AK4",  "PGK1",  "P4HA1",  "IGFBP5",  "ADM",  "PDK1",  "GSTM3",  "ALDOC",  "PFKFB4",  "GAPDH",  "BHLHE40",  "ANGPTL4",  "ENO1",  "STXBP1",  "SLC2A1",  "SPP1",  "FN1",  "DDIT4",  "MIF",  "C4orf3",  "TIMP3",  "HILPDA",  "ERO1A",  "PKM",  "MDK",  "SOX4",  "LDHA",  "EGLN1",  "SLC16A1",  "SLC2A3",  "SERPINB6",  "FBXO7",  "TPI1",  "ASPH",  "GPI",  "BNIP3L",  "PLOD2",  "RPS28",  "RNASE1",  "PLOD1",  "SLC7A5",  "XIST",  "TRIB2",  "FAM162A",  "NDRG1",  "WSB1",  "ANKRD37",  "KDM3A",  "RNF24",  "PRR34-AS1",  "BNIP3",  "EFR3B",  "ST13",  "RPL22",  "ZDBF2",  "HK2",  "FKBP9",  "SLC6A8",  "FXYD3",  "CFDP1",  "LEPROTL1",  "SRRM1",  "PLCG2",  "WDR45B",  "RPL36",  "BHLHE41",  "ALDOA",  "AMD1",  "RLF",  "DDX18",  "ZMYND8",  "HIST1H4C",  "NOL3",  "CITED1",  "A2M",  "UPP1",  "NQO1",  "TUT4",  "FXR1",  "TRIM2",  "DDX41",  "SARAF",  "ZRANB2",  "BRD2",  "FXYD1",  "HIST1H2AC",  "TF",  "TMSB10",  "IRF2BP2",  "SLC5A3",  "TNRC6B",  "INTU",  "P4HB",  "RNF19A")
geneSets <- GeneSet(genes, setName="Stress_hypoxia")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_malignant_subset_BT@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Stress_hypoxia<-getAUC(cells_AUC)
Stress_hypoxia<-t(Stress_hypoxia)
GC_malignant_subset_BT@meta.data<-cbind(GC_malignant_subset_BT@meta.data, Stress_hypoxia)
##########################        p_spec_B        ##########################
genes<-c("PSCA",  "TM4SF18",  "SLAMF1",  "ASB4",  "LINC02303",  "VAV3",  "GPRC5C",  "NTNG1",  "S100A7",  "PLCH1",  "HEPH",  "ZFP57",  "HHATL",  "SLC16A7",  "GHR",  "PSORS1C1",  "PRDM7",  "CSF2RA",  "ZDHHC2",  "CARD16",  "CA8",  "TLR1",  "KIF6",  "CYP39A1",  "PDE1A",  "ALDH1L2",  "PTCHD4",  "UGT2B7",  "GUCY1A2",  "SLITRK5",  "ITGA7",  "SYTL2",  "ENPP1",  "GPRC5A",  "CGNL1",  "CES1",  "SLC38A1",  "IL33",  "UNC5CL",  "RDH5",  "C5orf38",  "CDH1",  "FAM69A",  "SULT1A1",  "LINC02552",  "CYB5R2",  "RNF168",  "LINC01474",  "MET",  "FYB1",  "FABP7",  "POC1B",  "NNMT",  "LINC01315",  "PCLO",  "FKBP5",  "CPNE5",  "AC096564.1",  "AZGP1",  "PMP2",  "OAZ3",  "TDRD3",  "FAM27C",  "C15orf48",  "GXYLT2",  "AC090204.1",  "CXADR",  "CFI",  "RASSF2",  "CCDC71L",  "KCNN2",  "DPYD",  "PAWR",  "MIR29B2CHG",  "TRIM51",  "ASRGL1",  "HRK",  "AMOT",  "BCL2A1",  "CLIC2",  "SMPDL3A",  "HSPA4L",  "PERP",  "TBC1D7",  "BHLHE41",  "TMEM251",  "MRPL22",  "QPCT",  "ARL15",  "EIF3L",  "ARL6IP5",  "SUB1",  "LYPLAL1",  "SDHD",  "NDUFB4",  "OARD1",  "RPL22L1",  "FAM129A",  "TMEM14B",  "CCDC126")
geneSets <- GeneSet(genes, setName="p_spec_B")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_malignant_subset_BT@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
p_spec_B<-getAUC(cells_AUC)
p_spec_B<-t(p_spec_B)
GC_malignant_subset_BT@meta.data<-cbind(GC_malignant_subset_BT@meta.data, p_spec_B)
##########################        p_spec_A       ##########################
genes<-c("HAPLN1",  "WFDC2",  "FCER1G",  "AGT",  "NMU",  "IGDCC3",  "FBLN2",  "SEMA3F",  "MXRA5",  "NOS3",  "MYOC",  "FMN2",  "ARHGAP22",  "AMIGO2",  "OTULINL",  "TESC",  "ADRB1",  "ERVK-28",  "SHROOM3",  "LINC01597",  "CNDP1",  "STC2",  "DENND2A",  "PPP1R14A",  "EXOC3L2",  "C1QTNF4",  "CHST15",  "KCNJ12",  "DSC2",  "EMID1",  "SHISA2",  "SALL1",  "LINC02434",  "VSIR",  "TTYH1",  "ADGRF5",  "TNFRSF4",  "KLHL6",  "LRRC4B",  "TRIM71",  "FBLN1",  "TP53I11",  "RHOU",  "COL9A1",  "NKX2-5",  "NHSL2",  "PAPSS2",  "ATOH8",  "CCDC88C",  "TDRD12",  "C1QTNF3",  "CMTM8",  "SLC38A1",  "KIAA1217",  "FNDC5",  "IGF2BP1",  "DSG2",  "CBLN1",  "ABLIM1",  "SULT1C2",  "CRLF1",  "CIB2",  "PTGES3L",  "SULT1A2",  "LAMC3",  "TGM2",  "SULF2",  "LIN28B",  "CEACAM1",  "PXDN",  "DOK1",  "IMPA2",  "DUSP5",  "PLA2G16",  "TST",  "SERPINA5",  "AC012358.3",  "C1orf115",  "YJEFN3",  "DOK5",  "CLCN4",  "RFLNB",  "LINC00937",  "LINC00839",  "AC012651.1",  "NXN",  "KIAA1671",  "RASL10B",  "PDK4",  "HOMER3",  "ARSE",  "SULT1A1",  "ATP6V0A4",  "ITGA4",  "SLC19A3",  "Sep-06",  "DUSP15",  "PCLO",  "SDK2",  "CD68")
geneSets <- GeneSet(genes, setName="p_spec_A")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_malignant_subset_BT@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
p_spec_A<-getAUC(cells_AUC)
p_spec_A<-t(p_spec_A)
GC_malignant_subset_BT@meta.data<-cbind(GC_malignant_subset_BT@meta.data, p_spec_A)
VlnPlot(GC_malignant_subset_BT, "Mesenchymal")
FeatureScatter(GC_malignant_subset_BT, "SERPINB9", "Neural_crest_like")
ourscorelist <- c("Antigen_presenting","Interferon_responsive",  "Melanocytic", "Mesenchymal", "Mitochondrial", "Mitotic", "Neural_crest_like","p_spec_A","p_spec_B","Stress_hypoxia",  "Stress_paraspeckles", "MITF", "SOX10", "MLANA")

for (i in ourscorelist){
  p1 <- FeaturePlot(GC_malignant_subset_BT, features = i, cols = c("white", "black")) + NoAxes()
  pdf(paste0("/OUR_scores_UMAP_", i, ".pdf"), width = 5, height = 5)
  print(p1)
  dev.off()
  
}

####################################heatmap#########################################################
reglist <- c("Malignant_clusters", "Antigen_presenting","Interferon_responsive",  "Melanocytic", "Mesenchymal", "Mitochondrial", "Mitotic", "Neural_crest_like","p_spec_A","p_spec_B","Stress_hypoxia","Stress_paraspeckles")
matr <- FetchData(GC_malignant_subset_BT, vars = reglist, cells = NULL) # %>%
#apply(2, function(x) (x - min(x)) / (max(x) - min(x))) # %>%
#t
library(scales)
matr <- matr[ , !duplicated(colnames(matr))]
matr <- matr %>% group_by(Malignant_clusters) %>% summarise_all(mean)
#
Malignant_clusters <- levels(GC_malignant_subset_BT$Malignant_clusters)
Malignant_clusters
color_1 <- setNames(hue_pal()(11), Malignant_clusters)
color_list <- list(color_1 = Malignant_clusters) 
color_list
#names(color_list) <- list(Malignant_clusters_old)
Malignant_clusters <- matr$Malignant_clusters
ha = ComplexHeatmap::HeatmapAnnotation(Malignant_clusters = Malignant_clusters,
                                       annotation_name_side = "left",
                                       annotation_name_rot = 180,
                                       show_annotation_name = FALSE
                                       #col = color_list
)

#head(matr)
#matr <-apply(matr, 1, as.numeric)
head(matr)

matr$Malignant_clusters<- NULL
mat_num <- matrix(as.numeric(as.character((matr),    # Convert to numeric matrix
                                          ncol = ncol(matr))))
matr <- t(matr)
mat_scaled = t(apply(matr, 1, scale))
ht<-ComplexHeatmap::Heatmap(mat_scaled, 
                            cluster_rows = F, 
                            cluster_columns = FALSE,
                            show_column_names = FALSE,
                            top_annotation = ha,
                            #col = color_list,
                            column_split = Malignant_clusters)
pdf("/OUR_scores_heatmap_av.pdf", width = 10, height = 6)
draw(ht)
dev.off()



##########################        TOP_50_MES        ##########################
genes<-c("PDGFRB",  "NRP1",  "COL5A1",  "COL1A1",  "COL3A1",  "LUM",  "THY1",  "TCF4",  "ID3",  "PMEPA1",  "MRC2",  "DCN",  "COL1A2",  "COL4A1",  "MXRA8",  "TMSB4X",  "COL6A2",  "COL6A3",  "COL5A2",  "COL4A2",  "COL6A1",  "AQP1",  "FBLN1",  "IGFBP7",  "COL12A1",  "BGN",  "PCOLCE",  "VCAN",  "ANTXR2",  "PBX1",  "HTRA1",  "C1R",  "SELENOM",  "FBN1",  "SPARC",  "MARCKS",  "FN1",  "COL15A1",  "PALLD",  "CALD1",  "ITGA1",  "POSTN",  "CFH",  "PTN",  "ITGB1",  "H3F3B",  "TPM2",  "TIMP1",  "TMEM47",  "PRRX1")
geneSets <- GeneSet(genes, setName="TOP_50_MES")
geneSets
cells_rankings <- AUCell_buildRankings(GC_malignant_subset_BT@assays[["SCT"]]@counts,splitByBlocks=TRUE)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
TOP_50_MES<-getAUC(cells_AUC)
TOP_50_MES<-t(TOP_50_MES)
GC_malignant_subset_BT@meta.data<-cbind(GC_malignant_subset_BT@meta.data, TOP_50_MES)
pdf("/TOP_50_MES.pdf", width = 7.08, height = 5.8)
FeaturePlot(GC_malignant_subset_BT, features = "TOP_50_MES",label = T)
dev.off()
palette.malignant <- c("#77a600","#08519c","#ba5e45","#836ab1","grey", "#ff7f00" ,"#b30024", "#fcbd1f",  "#f28ac6", "#20b2aa","#e02887")
pdf("/TOP_50_MES_vln.pdf", width = 7, height = 6)
VlnPlot(GC_malignant_subset_BT, features = "TOP_50_MES", pt.size = 0.0, group.by = "Malignant_clusters", cols = palette.malignant) + theme(plot.margin=unit(c(1,1,1.5,1.5),"cm"), legend.position = 'none')
dev.off()
pdf("/Genes_across_patients_vln.pdf", width = 7, height = 6)
VlnPlot(GC_malignant_subset_BT, features = c("SOX10", "MITF", "CDH19", "S100A1", "TCF4"), pt.size = 0.001, group.by = "GC number") + theme(plot.margin=unit(c(1,1,1.5,1.5),"cm"), legend.position = 'none')
dev.off()
pdf("/TCF4_vln.pdf", width = 6, height = 6)
VlnPlot(GC_malignant_subset_BT, features = "TCF4", pt.size = 0, group.by = "Malignant_clusters") + theme(plot.margin=unit(c(1,1,1.5,1.5),"cm"), legend.position = 'none')
dev.off()

#saveRDS(GC_malignant_subset_BT, "/Malignant_BT.rds")
GC_malignant_subset_BT <- readRDS("/Malignant_BT.rds")

GC_malignant_subset_BT@meta.data$"Malignant_clusters" <- plyr::revalue(as.character(GC_malignant_subset_BT$seurat_clusters),
                                                                       c( "0" = "Melanocytic",
                                                                          "1" = "Mitochondrial(low_quality)",
                                                                          "2" = "Melanocytic",
                                                                          "3" = "Antigen_presentation",
                                                                          "4" = "Interferon_alpha_beta_response",
                                                                          "5" ="Stress (p53 response)",
                                                                          "6" ="Neural_Crest_like",
                                                                          "7" ="Stress (hypoxia response)",
                                                                          "8" ="Mitotic",
                                                                          "9" ="Patient_specific_A",
                                                                          "10" ="Mesenchymal_like",
                                                                          "11" ="Patient_specific_B"
                                                                       ))
############################# Alluvial Plots

cells <- rownames(GC_malignant_subset_BT@meta.data)
cells <- as.vector(cells)
library(data.table)
cnv_honey <- fread("/plotData_Honey_BADGER_for_score_immune.txt")
cnv_honey[1:5, 1:5]
dim(cnv_honey)
class(cnv_honey)
cnv_honey$V1 <- NULL
cnv_honey <- as.data.frame(cnv_honey)
cnv_honey[1:5, 1:5]

#subset cnv for only malignant cells
CNV_M <- dplyr::select(cnv_honey, all_of(cells))
dim(CNV_M)

CNV_M[1:5, 1:5]
class(CNV_M)
Malignant_clusters <- GC_malignant_subset_BT$Malignant_clusters
sample <- GC_malignant_subset_BT$orig.ident

ha = HeatmapAnnotation(Malignant_clusters = Malignant_clusters,
                       sample = sample,
                       annotation_name_side = "left",
                       annotation_name_rot = 00,
                       show_annotation_name = FALSE
)

mat_scaled = t(apply(CNV_M, 1, scale))
#col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
pdf("/CNV_heatmap.pdf", width = 10, height  = 10)
ht<-Heatmap(mat_scaled, 
            cluster_rows = TRUE, 
            cluster_columns = TRUE,
            show_column_names = FALSE,
            top_annotation = ha,
            clustering_method_rows = "ward.D",
            #col = col_fun,
            column_split = Malignant_clusters)
draw(ht)
dev.off()



#cluster the data

hc <- hclust(dist(t(CNV_M)), method='ward.D')
pdf("/dendrogram.pdf", width = 10, height  = 10)
plot(hc)
dev.off()
#prepare the output to add to seurat meta data
test <- as.data.frame(cutree(hc, k = 10)) # choose from the dendrogram K
test$CNV_clusters <- test$`cutree(hc, k = 10)`
test$`cutree(hc, k = 10)`<-NULL
test$CNV_clusters <-as.factor(test$CNV_clusters)
class(test$CNV_clusters)
class(test)
Idents(GC_malignant_subset_BT) <- "Malignant_clusters"

GC_malignant_subset_BT@meta.data$x_merge <- rownames(GC_malignant_subset_BT@meta.data)
test$x_merge <- rownames(test)
GC_malignant_subset_BT@meta.data <- GC_malignant_subset_BT@meta.data %>% inner_join(test, by="x_merge")
rownames(GC_malignant_subset_BT@meta.data) <- GC_malignant_subset_BT@meta.data$x_merge
#GC_malignant_subset_BT <- AddMetaData(GC_malignant_subset_BT, test, col.name="CNV_clusters")
#saveRDS(GC_malignant_subset_BT, "/Malignant_BT_cnv_clust.rds")
GC_malignant_subset_BT <- readRDS("/Malignant_BT_cnv_clust.rds")

### PLOT ####
custom_colors <- list()
colourCount = 11
getPalette = colorRampPalette(brewer.pal(11, "Dark2"))
Malignant_clusters <- getPalette(colourCount)



custom_colors$discrete <- c(Malignant_clusters)

# get sample and cluster names
GC_malignant_subset_BT@meta.data$Malignant_clusters <- as.factor(GC_malignant_subset_BT@meta.data$Malignant_clusters)
CNV_clusters <- levels(GC_malignant_subset_BT@meta.data$CNV_clusters)
Malignant_clusters <- levels(GC_malignant_subset_BT@meta.data$Malignant_clusters)

# create named vector holding the color assignments for both samples and
# Malignant_clusters
color_assignments <- setNames(
  c(custom_colors$discrete[1:length(CNV_clusters)], custom_colors$discrete[1:length(Malignant_clusters)]),
  c(CNV_clusters,Malignant_clusters)
)

data <- GC_malignant_subset_BT@meta.data %>%
  group_by(CNV_clusters, Malignant_clusters) %>%
  tally() %>%
  ungroup() %>%
  gather_set_data(1:2) %>%
  dplyr::mutate(
    x = factor(x, levels = unique(x)),
    y = factor(y, levels = unique(y))
  )

DataFrame(data)

data_labels <- tibble(
  group = c(
    rep('CNV_clusters', length(CNV_clusters)),
    rep('Malignant_clusters', length(Malignant_clusters))
  )
) %>%
  mutate(
    hjust = ifelse(group == 'CNV_clusters', 1, 0),
    nudge_x = ifelse(group == 'CNV_clusters', -0.1, 0.1)
  )

DataFrame(data_labels)

# create plot
p1 <- ggplot(data, aes(x, id = id, split = y, value = n)) +
  geom_parallel_sets(aes(fill = CNV_clusters), alpha = 0.75, axis.width = 0.15) +
  geom_parallel_sets_axes(aes(fill = y), color = 'black', axis.width = 0.1) +
  geom_text(
    aes(y = n, split = y), stat = 'parallel_sets_axes', fontface = 'bold',
    hjust = data_labels$hjust, nudge_x = data_labels$nudge_x
  ) +
  scale_x_discrete(labels = c('CNV_clusters','Malignant_clusters')) +
  scale_fill_manual(values = color_assignments) +
  theme_bw() +
  theme(
    legend.position = 'none',
    axis.title = element_blank(),
    axis.text.x = element_text(face = 'bold', colour = 'black', size = 15),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  )
pdf("/Alluvial_plot.pdf",  width = 9, height = 5.8)

p1
dev.off()


################ Signatures from Karras et al. paper ###########################################
##########################        Karras_Neural_Crest_like        ##########################
genes<-c("Tfap2b","Prnp","Mef2c","Gfra1","Cd200","Syt11","Thsd7a","Cxxc4","Sema5a","Tbx3","Kif26b","Efhd1","Neto2","Hmcn1","Igsf10","Olfml3","Rgs2","Nt5e","Morc4","Aqp1","Gfra2","Wnt4","Abcg2","Elovl5","Emilin1","Fibin")
genes <- convert_mouse_to_human_symbols(genes)
genes <- as.vector(genes)
genes
geneSets <- GeneSet(genes, setName="Karras_Neural_Crest_like")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_malignant_subset_BT@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
pdf("/Karras_Neural_Crest_like.pdf", width = 4, height = 4)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Karras_Neural_Crest_like<-getAUC(cells_AUC)
Karras_Neural_Crest_like<-t(Karras_Neural_Crest_like)
GC_malignant_subset_BT@meta.data<-cbind(GC_malignant_subset_BT@meta.data, Karras_Neural_Crest_like)
FeaturePlot(GC_malignant_subset_BT, features = "Karras_Neural_Crest_like", label = T) +NoAxes()
VlnPlot(GC_malignant_subset_BT, features = "Karras_Neural_Crest_like", pt.size = 0) + theme(legend.position = 'none')
VlnPlot(GC_malignant_subset_BT, features = "Karras_Neural_Crest_like", pt.size = 0) + theme(legend.position = 'none')
VlnPlot(GC_malignant_subset_BT, features = "Karras_Neural_Crest_like", pt.size = 0, group.by = "orig.ident") + theme(legend.position = 'none')
VlnPlot(GC_malignant_subset_BT, features = "Karras_Neural_Crest_like", pt.size = 0, group.by = "GC number") + theme(legend.position = 'none')
dev.off()

##########################        Karras_Mesenchymal_like        ##########################
genes<-c("Fap","Lama2","Prrx1","Slit3","Abi3bp","Loxl1","Cdh11","Col6a3","Col6a2","Loxl2","Mfap5","Fbn1","Col4a2","Pcolce","Lum","Col4a1","Col5a2","Thy1","Fbn2","Bgn","Tgfbi","Pdgfrb","Sulf1","Inhba","Col1a1","Col3a1","Col1a2","Ctsk","Dcn","Serpinf1","Sparc","Fstl1")
genes <- convert_mouse_to_human_symbols(genes)
genes <- as.vector(genes)
genes
geneSets <- GeneSet(genes, setName="Karras_Mesenchymal_like")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_malignant_subset_BT@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
pdf("/Karras_Mesenchymal_like.pdf", width = 7.08, height = 5.8)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Karras_Mesenchymal_like<-getAUC(cells_AUC)
Karras_Mesenchymal_like<-t(Karras_Mesenchymal_like)
GC_malignant_subset_BT@meta.data<-cbind(GC_malignant_subset_BT@meta.data, Karras_Mesenchymal_like)
FeaturePlot(GC_malignant_subset_BT, features = "Karras_Mesenchymal_like", label = T)
VlnPlot(GC_malignant_subset_BT, features = "Karras_Mesenchymal_like", pt.size = 0) + theme(legend.position = 'none')
VlnPlot(GC_malignant_subset_BT, features = "Karras_Mesenchymal_like", pt.size = 0) + theme(legend.position = 'none')
VlnPlot(GC_malignant_subset_BT, features = "Karras_Mesenchymal_like", pt.size = 0, group.by = "orig.ident") + theme(legend.position = 'none')
VlnPlot(GC_malignant_subset_BT, features = "Karras_Mesenchymal_like", pt.size = 0, group.by = "GC number") + theme(legend.position = 'none')
dev.off()

##########################        Karras_Melanocytic        ##########################
genes<-c("Mlph","Mitf","Mlana","Pmel","Slc45a2","Apoe","Dct","Gpnmb","Tyr","Cited1","Bace2","Cox17","Cox7a2","Ndufb2","Ndufb4","Uqcr10","Uqcrb")
genes <- convert_mouse_to_human_symbols(genes)
genes <- as.vector(genes)
genes
geneSets <- GeneSet(genes, setName="Karras_Melanocytic")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_malignant_subset_BT@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
pdf("/Karras_Melanocytic.pdf", width = 7.08, height = 5.8)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Karras_Melanocytic<-getAUC(cells_AUC)
Karras_Melanocytic<-t(Karras_Melanocytic)
GC_malignant_subset_BT@meta.data<-cbind(GC_malignant_subset_BT@meta.data, Karras_Melanocytic)
FeaturePlot(GC_malignant_subset_BT, features = "Karras_Melanocytic", label = T)
VlnPlot(GC_malignant_subset_BT, features = "Karras_Melanocytic", pt.size = 0) + theme(legend.position = 'none')
VlnPlot(GC_malignant_subset_BT, features = "Karras_Melanocytic", pt.size = 0) + theme(legend.position = 'none')
VlnPlot(GC_malignant_subset_BT, features = "Karras_Melanocytic", pt.size = 0, group.by = "orig.ident") + theme(legend.position = 'none')
VlnPlot(GC_malignant_subset_BT, features = "Karras_Melanocytic", pt.size = 0, group.by = "GC number") + theme(legend.position = 'none')
dev.off()



##########################        Karras_RNA_processing        ##########################
genes<-c("Pabpn1","Ccnl1","Pnisr","Sfpq","Ccnl2","Bclaf1","Hnrnpu","Rbm25","Luc7l2","Ddx17","Cpsf6","Snrnp70","Srek1","Prpf4b","Srrm2","Rbm39","Son","Ddx5","Prpf38b","Pan3","Tra2a","Hnrnph1","Fus")
genes <- convert_mouse_to_human_symbols(genes)
genes <- as.vector(genes)
genes
geneSets <- GeneSet(genes, setName="Karras_RNA_processing")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_malignant_subset_BT@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
pdf("/Karras_RNA_processing.pdf", width = 7.08, height = 5.8)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Karras_RNA_processing<-getAUC(cells_AUC)
Karras_RNA_processing<-t(Karras_RNA_processing)
GC_malignant_subset_BT@meta.data<-cbind(GC_malignant_subset_BT@meta.data, Karras_RNA_processing)
FeaturePlot(GC_malignant_subset_BT, features = "Karras_RNA_processing", label = T)
VlnPlot(GC_malignant_subset_BT, features = "Karras_RNA_processing", pt.size = 0) + theme(legend.position = 'none')
VlnPlot(GC_malignant_subset_BT, features = "Karras_RNA_processing", pt.size = 0) + theme(legend.position = 'none')
VlnPlot(GC_malignant_subset_BT, features = "Karras_RNA_processing", pt.size = 0, group.by = "orig.ident") + theme(legend.position = 'none')
VlnPlot(GC_malignant_subset_BT, features = "Karras_RNA_processing", pt.size = 0, group.by = "GC number") + theme(legend.position = 'none')
dev.off()


##########################        Karras_Stress_hypoxia        ##########################
genes<-c("Bnip3","Tpi1","Slc2a1","Mif","Vldlr","Hk2","Vegfa","Ldha","Pfkl","P4ha1","Fam162a","Bhlhe40","Pgk1","Aldoa","Pfkp","Pdk1","Hspa9","Pgam1","Pkm","Atf4")
genes <- convert_mouse_to_human_symbols(genes)
genes <- as.vector(genes)
genes
geneSets <- GeneSet(genes, setName="Karras_Stress_hypoxia")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_malignant_subset_BT@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
pdf("/Karras_Stress_hypoxiaR.pdf", width = 7.08, height = 5.8)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Karras_Stress_hypoxia<-getAUC(cells_AUC)
Karras_Stress_hypoxia<-t(Karras_Stress_hypoxia)
GC_malignant_subset_BT@meta.data<-cbind(GC_malignant_subset_BT@meta.data, Karras_Stress_hypoxia)
FeaturePlot(GC_malignant_subset_BT, features = "Karras_Stress_hypoxia", label = T)
VlnPlot(GC_malignant_subset_BT, features = "Karras_Stress_hypoxia", pt.size = 0) + theme(legend.position = 'none')
VlnPlot(GC_malignant_subset_BT, features = "Karras_Stress_hypoxia", pt.size = 0) + theme(legend.position = 'none')
VlnPlot(GC_malignant_subset_BT, features = "Karras_Stress_hypoxia", pt.size = 0, group.by = "orig.ident") + theme(legend.position = 'none')
VlnPlot(GC_malignant_subset_BT, features = "Karras_Stress_hypoxia", pt.size = 0, group.by = "GC number") + theme(legend.position = 'none')
dev.off()


##########################        Karras_Stem_like         ##########################
genes<-c("Gtse1","Cep170b","Tap1","Ercc5","Kank3","Rap2a","Ei24","Zmat3","Fas","Igdcc4","Notch3","Vcan","Mybl1","Dusp15","Dgka","Nrp1","Slc27a3","Slc19a2","Dgki","Siva1","Rbm38","Cad","Mcm7","Nme4","Thyn1","Gpx3","Arl4c","Ass1","Efnb2","Frrs1","Sdc1","Itga6","Icam1","Chst3","Nes","Sox4","Fat1","Angptl2")
genes <- convert_mouse_to_human_symbols(genes)
genes <- as.vector(genes)
genes
geneSets <- GeneSet(genes, setName="Karras_Stem_like")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_malignant_subset_BT@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
pdf("/Karras_Stem_like.pdf", width = 4, height = 4)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Karras_Stem_like<-getAUC(cells_AUC)
Karras_Stem_like<-t(Karras_Stem_like)
GC_malignant_subset_BT@meta.data<-cbind(GC_malignant_subset_BT@meta.data, Karras_Stem_like)
FeaturePlot(GC_malignant_subset_BT, features = "Karras_Stem_like", label = T) +NoAxes()
VlnPlot(GC_malignant_subset_BT, features = "Karras_Stem_like", pt.size = 0) + theme(legend.position = 'none')
VlnPlot(GC_malignant_subset_BT, features = "Karras_Stem_like", pt.size = 0, group.by = "orig.ident") + theme(legend.position = 'none')
VlnPlot(GC_malignant_subset_BT, features = "Karras_Stem_like", pt.size = 0, group.by = "GC number") + theme(legend.position = 'none')
dev.off()



##########################        Karras_Antigen_presenting        ##########################
genes<-c("Gbp4","Gbp6","Irf1","Nlrc5","Stat2","Cd74","Il12rb1","B2m","Gbp2","Ifit3","Gbp3","Psmb8","Stat1","Ifit1","Isg15","Psmb9","Gbp7","Xaf1","Tap1","Gbp5","Usp18","Ciita","Irf7","Ube2l6","Tap2","Psme1","Ifitm3","Psme2","Tapbp","Ifi35","Eif2ak2","Irf9","Ddx58","Trim25","Ifit2","Irf8","H2-DMa","H2-Aa","H2-DMb1","H2-Eb1","H2-Ab1","H2-T23")
genes <- convert_mouse_to_human_symbols(genes)
genes <- as.vector(genes)
genes
 genes <-unique(genes)
 genes
geneSets <- GeneSet(genes, setName="Karras_Antigen_presenting")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_malignant_subset_BT@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
pdf("/Karras_Antigen_presenting.pdf", width = 7.08, height = 5.8)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Karras_Antigen_presenting<-getAUC(cells_AUC)
Karras_Antigen_presenting<-t(Karras_Antigen_presenting)
GC_malignant_subset_BT@meta.data<-cbind(GC_malignant_subset_BT@meta.data, Karras_Antigen_presenting)
FeaturePlot(GC_malignant_subset_BT, features = "Karras_Antigen_presenting", label = T)
VlnPlot(GC_malignant_subset_BT, features = "Karras_Antigen_presenting", pt.size = 0) + theme(legend.position = 'none')
VlnPlot(GC_malignant_subset_BT, features = "Karras_Antigen_presenting", pt.size = 0, group.by = "orig.ident") + theme(legend.position = 'none')
VlnPlot(GC_malignant_subset_BT, features = "Karras_Antigen_presenting", pt.size = 0, group.by = "GC number") + theme(legend.position = 'none')
dev.off()
VlnPlot(GC_malignant_subset_BT, features = c("NGFR", "NRXN", "SNCA", "SYT11", "GPHN"), pt.size = 0, group.by = "Malignant_clusters") + theme(legend.position = 'none')


############## heatmap 
####################################heatmap#########################################################
reglist <- c("Malignant_clusters","Karras_Antigen_presenting","Karras_Melanocytic", "Karras_Mesenchymal_like", "Karras_Neural_Crest_like","Karras_Stress_hypoxia","Karras_RNA_processing","Karras_Stem_like")
matr <- FetchData(GC_malignant_subset_BT, vars = reglist, cells = NULL) # %>%
#apply(2, function(x) (x - min(x)) / (max(x) - min(x))) # %>%
#t
library(scales)
matr <- matr[ , !duplicated(colnames(matr))]
matr <- matr %>% group_by(Malignant_clusters) %>% summarise_all(mean)
#
Malignant_clusters <- levels(GC_malignant_subset_BT$Malignant_clusters)
Malignant_clusters
color_1 <- setNames(hue_pal()(11), Malignant_clusters)
color_list <- list(color_1 = Malignant_clusters) 
color_list
#names(color_list) <- list(Malignant_clusters_old)
Malignant_clusters <- matr$Malignant_clusters
ha = ComplexHeatmap::HeatmapAnnotation(Malignant_clusters = Malignant_clusters,
                                       annotation_name_side = "left",
                                       annotation_name_rot = 180,
                                       show_annotation_name = FALSE
                                      # col = color_list
)

#head(matr)
#matr <-apply(matr, 1, as.numeric)
head(matr)

matr$Malignant_clusters <- NULL
mat_num <- matrix(as.numeric(as.character((matr),    # Convert to numeric matrix
                                          ncol = ncol(matr))))
matr <- t(matr)
mat_scaled = t(apply(matr, 1, scale)) #by row
ht<-ComplexHeatmap::Heatmap(mat_scaled, 
                            cluster_rows = FALSE, 
                            cluster_columns = FALSE,
                            show_column_names = FALSE,
                            top_annotation = ha,
                            #col = color_list,
                            column_split = Malignant_clusters)
pdf("/scores_heatmap_av_Karras_et_al.pdf", width = 10, height = 6)
draw(ht)
dev.off()









