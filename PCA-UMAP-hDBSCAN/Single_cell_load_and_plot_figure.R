library("Seurat")
#library(biomaRt)
library(dplyr)
library(tidyr)
library("ggplot2")
library(randomcoloR)
library(RColorBrewer)
library("monocle")
library(dbscan)

###########please download files from figshare: https://figshare.com/s/d3e33d8a53163737481d


setwd("CRCLM-master/PCA-UMAP-hDBSCAN")
###########load files
pbmc<-readRDS("20210228_14sample_fil_Bcell_pbmc_reclustering.rds") #this file can be find in figureshare or Zenodo
pbmc_pure<-readRDS("20210228_14sample_fil_Bcell_pbmc_reclusteringpbmc_pure.rds") #this file can be find in figureshare or Zenodo
pbmc_pure_1<-readRDS("20210228_14sample_fil_Bcell_pbmc_reclusteringpbmc_pure_1.rds") #this file can be find in figureshare or Zenodo
pbmc_part_metastatic_R<-readRDS("20210228_14sample_fil_Bcell_pbmc_reclusteringpbmc_pure_1_half.rds") #this file can be find in figureshare or Zenodo
load("Pseudotime_monocle_objB_cell1") 

###########plot
palette_seurat_clusters <- randomColor(length(unique(pbmc$seurat_clusters)))
palette_LR_MN <- brewer.pal( length(unique(pbmc$LR_MN)),name ="Set3")
palette_sample <- randomColor(length(unique(pbmc$sample)))
palette_meta<- randomColor(length(unique(pbmc$meta)))
palette_seurat_clusters <- randomColor(length(unique(pbmc$seurat_clusters)))

pdf(file=paste0("Figure_",Sys.Date(),"B_cell_DimPlot_dbscan.pdf"))
print(DimPlot(pbmc, reduction = "umap",cols=palette_seurat_clusters,pt.size=0.3))
dev.off()
pdf(file=paste0("Figure_",Sys.Date(),"B_cell_DimPlot_dbscan_puro.pdf"))
palette_seurat_clusters <- randomColor(length(unique(pbmc_pure$seurat_clusters)))
print(DimPlot(pbmc_pure, reduction = "umap",cols=palette_seurat_clusters,pt.size=0.3))
dev.off()
pdf(file=paste0("Figure_",Sys.Date(),"B_cell_DimPlot_meta.pdf"),6,5)
print(DimPlot(pbmc, reduction = "umap", cols=palette_meta, group.by = "meta",pt.size=0.3))
dev.off()


pbmc.markers <- FindAllMarkers(pbmc_pure, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, num.cores = 64)

#pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, num.cores = 64)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
pdf(file=paste0("Figure_",Sys.Date(),"B_cellS6_DoHeatmap_dbscan.pdf"),10,15)
print(DoHeatmap(pbmc_pure, features = top10$gene) + NoLegend())
dev.off()

#################plot Hdbscan clusters percentage figure 3 F

seurat_clusters_count<-data.frame(seurat_clusters=pbmc$seurat_clusters,Group=pbmc$meta)

result_race<-seurat_clusters_count %>% 
  group_by(Group) %>% count(seurat_clusters)


result_race_O<-result_race
result_race$n <- ave(result_race$n, result_race$seurat_clusters, FUN=function(x) x/max(x))
result_race$value <-result_race_O$n

###########


gender_plot<-ggplot(result_race, aes(reorder(seurat_clusters,n),fill=as.factor(Group), y=n)) + 
  geom_bar(position="fill", stat="identity")+ scale_fill_manual(values = palette_LR_MN)+  
  geom_text(aes(label=value),stat='identity',position=position_fill(vjust=0.5))+
  labs(fill = "mate",x="Hdbscan_clusters")+ ggtitle("B cell") +
  coord_flip()



pdf(file=paste0("Figure_",Sys.Date(),"B_cells_geom_bar_dbscaLR_mate.pdf"))
gender_plot
dev.off()

#################plot Hdbscan clusters percentage figure

seurat_clusters_count<-data.frame(seurat_clusters=pbmc$seurat_clusters,Group=pbmc$LR_MN)

result_race<-seurat_clusters_count %>% 
  group_by(Group) %>% count(seurat_clusters)


result_race_O<-result_race
result_race$n <- ave(result_race$n, result_race$seurat_clusters, FUN=function(x) x/max(x))
result_race$value <-result_race_O$n

###########


gender_plot<-ggplot(result_race, aes(reorder(seurat_clusters,n),fill=as.factor(Group), y=n)) + 
  geom_bar(position="fill", stat="identity")+ scale_fill_manual(values = palette_LR_MN)+  
  geom_text(aes(label=value),stat='identity',position=position_fill(vjust=0.5))+
  labs(fill = "LR_MN",x="Hdbscan_clusters")+ ggtitle("B cell") +
  coord_flip()



pdf(file=paste0("Figure_",Sys.Date(),"B_cells_geom_bar_dbscaLR_MNn.pdf"))
gender_plot
dev.off()


pdf(file=paste0("Figure_",Sys.Date(),"B_cell_DimPlot_dbscan_puro.pdf"))
palette_seurat_clusters <- randomColor(length(unique(pbmc_pure$seurat_clusters)))
print(DimPlot(pbmc_pure_1, reduction = "umap",cols=palette_seurat_clusters,pt.size=0.3))
dev.off()

pdf(file=paste0("Figure_",Sys.Date(),"B_cell_culs1_DimPlot_LR_MN.pdf"),7,5)
print(DimPlot(pbmc_pure_1, reduction = "umap", cols=palette_LR_MN, group.by = "LR_MN",pt.size=1))
dev.off()




seurat_clusters_count<-data.frame(seurat_clusters=pbmc_pure_1$seurat_clusters,Group=pbmc_pure_1$meta)

result_race<-seurat_clusters_count %>% 
  group_by(Group) %>% count(seurat_clusters)


result_race_O<-result_race

result_race$n <- ave(result_race$n, result_race$seurat_clusters, FUN=function(x) x/max(x)) 
#result_race$n <- ave(result_race$n, result_race$seurat_clusters, FUN=function(x) x/max(x))
result_race$value <-result_race_O$n

###########


gender_plot<-ggplot(result_race, aes(reorder(seurat_clusters,n),fill=as.factor(Group), y=n)) + 
  geom_bar(position="fill", stat="identity")+ scale_fill_manual(values = palette_LR_MN)+  
  geom_text(aes(label=value),stat='identity',position=position_fill(vjust=0.5))+
  labs(fill = "Group",x="Clusters")+ ggtitle("B cell") +
  coord_flip()



pdf(file=paste0("Figure_",Sys.Date(),"B_cells_geom_bar_kmeans_LR_MN.pdf"))
gender_plot
dev.off()


################### percentage 


seurat_clusters_count<-data.frame(seurat_clusters=pbmc_pure_1$seurat_clusters,Group=pbmc_pure_1$LR_MN)

result_race<-seurat_clusters_count %>% 
  group_by(Group) %>% count(seurat_clusters)


result_race_O<-result_race

result_race$n <- ave(result_race$n, result_race$seurat_clusters, FUN=function(x) x/max(x)) 
#result_race$n <- ave(result_race$n, result_race$seurat_clusters, FUN=function(x) x/max(x))
result_race$value <-result_race_O$n

###########


gender_plot<-ggplot(result_race, aes(reorder(seurat_clusters,n),fill=as.factor(Group), y=n)) + 
  geom_bar(position="fill", stat="identity")+ scale_fill_manual(values = palette_LR_MN)+  
  geom_text(aes(label=value),stat='identity',position=position_fill(vjust=0.5))+
  labs(fill = "Group",x="Clusters")+ ggtitle("B cell") +
  coord_flip()



pdf(file=paste0("FigureS3D_",Sys.Date(),"B_cells_geom_bar_kmeans_LR_MN.pdf"))
gender_plot
dev.off()



###########monocle analysis


patient_name<-"B_cell" 
Selected_population <- 1

gc()


pbmc_part_metastatic_R <- pbmc_pure_1[,which(pbmc_pure_1$seurat_clusters ==Selected_population)]

matrix_temp <- pbmc_part_metastatic_R@assays[["RNA"]]@counts
matrix_temp_scale <- pbmc_part_metastatic_R@assays[["RNA"]]@data

gc()

pdf(file=paste0("Figure_",Sys.Date(),"B_cell_plot_cell_trajectory_dbscan_puro1_meta.pdf"))
print(plot_cell_trajectory(HSMM_myo, color_by= "meta")+ scale_color_manual(breaks = waiver(),values=palette_meta))
dev.off()


plot_gene <- "IGHG4"
pdf(file=paste0("Figure_",Sys.Date(),"B_cell_plot_cell_trajectory_dbscan_puro1_",plot_gene,".pdf"))
plot_cell_trajectory(HSMM_myo, color_by =(matrix_temp_scale[plot_gene,]))+scale_colour_gradientn(colours = c("#2166ac","#67a9cf","#d1e5f0","#fddbc7","#ef8a62","#b2182b"),limits=c(min(matrix_temp_scale[plot_gene,]),max(matrix_temp_scale[plot_gene,])))+labs(title=paste(plot_gene))
dev.off()






###########differential experssion analysis

expressed_genes <- row.names(subset(fData(HSMM_myo), num_cells_expressed >=10))

diff_test_res <- differentialGeneTest(HSMM_myo[expressed_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)",cores=8)
sig_gene_names <- row.names(subset(diff_test_res, (qval < 5**-100 & use_for_ordering=="TRUE")))



print(plot_pseudotime_heatmap(HSMM_myo[sig_gene_names,],
                              num_clusters = 6,
                              cores = 8,
                              show_rownames = T))


marker_gene<-FindMarkers(pbmc_pure,1)


marker_gene <- intersect(rownames(HSMM_myo),rownames(marker_gene))


# figure
print(plot_pseudotime_heatmap(HSMM_myo[marker_gene,],
                              num_clusters = 5,
                              cores = 8,
                              show_rownames = T))



TF_genes <- c("NME2","EVI1","COMP1","CEBPB","SYNCRIP","GABP","TAXCREB","NPM1","CTR9","GTF2E2","ARHGAP35","YY1","TAF9B","HJURP","NEUROG3","SOX11","MAPK3","PSMB5","ZNF318","GTF2A2","FOXR2","PER1","ZZZ3")


TF_genes  <- intersect(rownames(HSMM_myo),TF_genes )



# figure
print(plot_pseudotime_heatmap(HSMM_myo[TF_genes,],
                              num_clusters = 3,
                              cores = 8,
                              show_rownames = T))



