library("Seurat")
library(biomaRt)
library(dplyr)
library(tidyr)
library("ggplot2")
library(randomcoloR)
library(RColorBrewer)
library("monocle")
library("ggpubr")


setwd("E:/XXXXXXXX/S_aggr12/S_aggr12/14sample_fil")

#load data
#######################################################
pbmc<-readRDS("14sample_fil_Bcell_pbmc_reclusteringpbmc_pure.rds")


#samples number
#######################################################
all_cells_total<-data.frame(Barcode=colnames(pbmc),(data.frame(x=colnames(pbmc)) %>% separate(x, c("A", "B"))))
sample<-recode(all_cells_total$B, "1"="S1","2"="S2","3"="S3","4"="S10","5"="S5","6"="S6","7"="S7","8"="S8","9"="S9","10"="S11","11"="S12","12"="S13","13"="S14")
all_cells_total<-data.frame(all_cells_total,sample)
meta<-recode(all_cells_total$sample, "S1"="no-metastasis","S2"="no-metastasis","S3"="liver-metastasis","S10"="liver-metastasis","S5"="liver-metastasis","S6"="no-metastasis","S7"="liver-metastasis","S8"="no-metastasis","S9"="liver-metastasis","S11"="no-metastasis","S12"="no-metastasis","S13"="liver-metastasis","S14"="no-metastasis")
side<-recode(all_cells_total$sample, "S1"="metastatic-R","S2"="metastatic-R","S3"="metastatic-R","S10"="metastatic-L","S5"="metastatic-R","S6"="metastatic-R","S7"="metastatic-R","S8"="metastatic-R","S9"="metastatic-L","S11"="metastatic-L","S12"="metastatic-L","S13"="metastatic-L","S14"="metastatic-L")
pbmc@meta.data[["sample"]] <-sample
pbmc@meta.data[["meta"]] <-meta
pbmc@meta.data[["side"]] <-side
pbmc@meta.data[["LR_MN"]]<-(paste(pbmc@meta.data[["side"]],pbmc@meta.data[["meta"]],sep = "_"))


#Seurat step (PCA-UMAP)
#######################################################
gc()
pbmc <- ScaleData(pbmc)
pbmc <- ScaleData(pbmc, vars.to.regress = c("nFeature_RNA","nCount_RNA", "percent.mt"), features = rownames(pbmc))
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
unique(pbmc$orig.ident)
pbmc_list <- SplitObject(pbmc, split.by = "orig.ident")
object.anchors <- FindIntegrationAnchors(object.list = pbmc_list, dims = 1:5)
pbmc <- IntegrateData(anchorset = object.anchors, dims = 1:5)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc)
pbmc <-RunUMAP(pbmc, dims = 1:20,min.dist = 0.0, n.neighbors = 200L)
pbmc <- RunTSNE(pbmc, dims = 1:10)
pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 0.3)


#Plot PCA-UMAP reslut
#######################################################
palette_seurat_clusters <- randomColor(length(unique(pbmc$seurat_clusters)))
palette_LR_MN <- brewer.pal( length(unique(pbmc$LR_MN)),name ="Set3")
palette_sample <- randomColor(length(unique(pbmc$sample)))
palette_meta<- randomColor(length(unique(pbmc$meta)))



pdf(file=paste0("Figure_",Sys.Date(),"B_cell_DimPlot_clusters.pdf"))
print(DimPlot(pbmc, reduction = "umap",cols=palette_seurat_clusters,pt.size=0.3))
dev.off()
pdf(file=paste0("Figure_",Sys.Date(),"B_cell_DimPlot_LR_MN.pdf"),7,5)
print(DimPlot(pbmc, reduction = "umap", cols=palette_LR_MN, group.by = "LR_MN",pt.size=0.3))
dev.off()
pdf(file=paste0("Figure_",Sys.Date(),"B_cell_DimPlot_sample.pdf"))
print(DimPlot(pbmc, reduction = "umap",cols=palette_sample, group.by = "sample",pt.size=0.3))
dev.off()
pdf(file=paste0("Figure_",Sys.Date(),"B_cell_DimPlot_meta.pdf"),6,5)
print(DimPlot(pbmc, reduction = "umap", cols=palette_meta, group.by = "meta",pt.size=0.3))
dev.off()



#hdbscan reclustering
#######################################################
library(dbscan)


umap_data <- data.frame(pbmc@reductions[["umap"]]@cell.embeddings)
dbscan_result <-hdbscan(pbmc@reductions[["umap"]]@cell.embeddings,minPts=80)
sortcluster<-sort(table(dbscan_result$cluster),decreasing = T)

orgi_num<-names(sortcluster)
target_num<-c(0:9)

long<- c()
for (i in seq_along(sortcluster))
{
  #i=1
  tep<-noquote(paste0(dQuote(orgi_num[i]), noquote("="),dQuote(target_num[i])))
  
  long<-paste(tep,long,sep=",")
}


recoded_dbscan<-recode(dbscan_result$cluster,"5"="8","4"="7","0"="6","7"="5","2"="4","1"="3","3"="2","6"="1","8"="0")

n <- length(unique(dbscan_result$cluster ))
palette <- distinctColorPalette(n)

ggplot(umap_data,mapping = aes(x= UMAP_1, y= UMAP_2, color= as.character(recoded_dbscan) ))+
  geom_point(size = 1)+
  scale_color_manual(values =  palette)

ident<-as.factor(recoded_dbscan)
names(ident) <-colnames(pbmc)

pbmc@active.ident <- ident
pbmc$seurat_clusters <- recoded_dbscan


#Filter poorly clustered clusters
#######################################################
pbmc_pure<-subset(pbmc, seurat_clusters !=6 )



#Plot PCA-UMAP-Hdbscan reslut
#######################################################
palette_seurat_clusters <- randomColor(length(unique(pbmc$seurat_clusters)))

pdf(file=paste0("Figure_",Sys.Date(),"B_cell_DimPlot_dbscan.pdf"))
print(DimPlot(pbmc, reduction = "umap",cols=palette_seurat_clusters,pt.size=0.3))
dev.off()
pdf(file=paste0("Figure_",Sys.Date(),"B_cell_DimPlot_dbscan_puro.pdf"))
palette_seurat_clusters <- randomColor(length(unique(pbmc_pure$seurat_clusters)))
print(DimPlot(pbmc_pure, reduction = "umap",cols=palette_seurat_clusters,pt.size=0.3))
dev.off()



#Plot genes on PCA-UMAP-Hdbscan reslut
#######################################################
#plot_gene <- "CD79A"

pdf(file=paste0("Figure_",Sys.Date(),"B_cell_FeaturePlot_puro1_",plot_gene,".pdf"))
print(FeaturePlot(pbmc_pure, features = plot_gene,cols = c("#2166ac","#67a9cf","#d1e5f0","#fddbc7","#ef8a62","#b2182b") ))
dev.off()


pdf(file=paste0("Figure_",Sys.Date(),"B_cell_VlnPlot_puro1_",plot_gene,".pdf"))
print(VlnPlot(object = pbmc_pure, features =plot_gene))
dev.off()


#Differential expression analysis
#######################################################
pbmc.markers <- FindAllMarkers(pbmc_pure, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, num.cores = 64)
#pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, num.cores = 64)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

pdf(file=paste0("Figure_",Sys.Date(),"B_cellS6_DoHeatmap_dbscan.pdf"),10,15)
print(DoHeatmap(pbmc_pure, features = top10$gene) + NoLegend())
dev.off()




#####################################################################################################


#load data
#######################################################
pbmc_pure_1<-readRDS("20210228_14sample_fil_Bcell_pbmc_reclusteringpbmc_pure_1.rds")


#Kmean clustering
#######################################################
umap_data <- data.frame(pbmc_pure_1@reductions[["umap"]]@cell.embeddings)

kmeans_result <-kmeans(pbmc_pure_1@reductions[["umap"]]@cell.embeddings,2)


pbmc_pure_1@active.ident <- as.factor(kmeans_result$cluster)
pbmc_pure_1$seurat_clusters <- as.factor(kmeans_result$cluster)

palette_seurat_clusters <- randomColor(length(unique(pbmc_pure_1$seurat_clusters)))

pdf(file=paste0("Figure_",Sys.Date(),"B_cell_Kmean_DimPlot_dbscan.pdf"))
print(DimPlot(pbmc_pure_1, reduction = "umap",cols=palette_seurat_clusters,pt.size=0.3))
dev.off()
pdf(file=paste0("Figure_",Sys.Date(),"B_cell_Kmean_DimPlot_dbscan_puro.pdf"))
palette_seurat_clusters <- randomColor(length(unique(pbmc_pure$seurat_clusters)))
print(DimPlot(pbmc_pure_1, reduction = "umap",cols=palette_seurat_clusters,pt.size=0.3))
dev.off()





#Monocle analysis
#######################################################
patient_name<-"B_cell" 
Selected_population <- 1

gc()
pbmc_part_metastatic_R <- pbmc_pure_1[,which(pbmc_pure_1$seurat_clusters ==Selected_population)]
matrix_temp <- pbmc_part_metastatic_R@assays[["RNA"]]@counts
matrix_temp_scale <- pbmc_part_metastatic_R@assays[["RNA"]]@data

gc()

HSMM <- newCellDataSet(as.matrix(matrix_temp))
HSMM@phenoData@data[["meta"]] <- pbmc_part_metastatic_R@meta.data[["meta"]]
HSMM@phenoData@data[["sample"]] <- pbmc_part_metastatic_R@meta.data[["sample"]]
HSMM@phenoData@data[["LR_MN"]] <- pbmc_part_metastatic_R@meta.data[["LR_MN"]]

HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)
HSMM <- detectGenes(HSMM, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 10))
HSMM <- detectGenes(HSMM, min_expr = 0.1)
L <- log(exprs(HSMM[expressed_genes,]))
melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))
gc()

disp_table <- dispersionTable(HSMM)
ordering_genes <- subset(disp_table,
                         mean_expression >= 0.01 &
                           dispersion_empirical >= 1 * dispersion_fit)$gene_id


HSMM <- setOrderingFilter(HSMM, ordering_genes)
gc()
HSMM <- reduceDimension(HSMM, max_components=2, num_dim = 20,
                        reduction_method = 'tSNE', verbose = T)
HSMM <- clusterCells(HSMM, num_clusters=2)
gc()
HSMM_myo <- estimateDispersions(HSMM)
disp_table <- dispersionTable(HSMM_myo)
ordering_genes <- subset(disp_table,
                         mean_expression >= 0.01 &
                           dispersion_empirical >=1 * dispersion_fit)$gene_id
HSMM_myo <- setOrderingFilter(HSMM_myo, ordering_genes)

gc()
HSMM_myo <- reduceDimension(HSMM_myo, max_components=2)
HSMM_myo <- orderCells(HSMM_myo)

#save(HSMM_myo,file=paste0("Pseudotime_monocle_obj",patient_name,Selected_population))

print(plot_cell_trajectory(HSMM_myo))
print(plot_cell_trajectory(HSMM_myo, color_by="Pseudotime"))
print(plot_cell_trajectory(HSMM_myo, color_by= "side")) 


pdf(file=paste0("Figure2B_",Sys.Date(),"B_cell_plot_cell_trajectory_dbscan_puro1_meta.pdf"))
print(plot_cell_trajectory(HSMM_myo, color_by= "meta")+ scale_color_manual(breaks = waiver(),values=palette_meta))
dev.off()

pdf(file=paste0("Figure2B_",Sys.Date(),"B_cell_plot_cell_trajectory_dbscan_puro1_LR_MN.pdf"))
print(plot_cell_trajectory(HSMM_myo, color_by= "LR_MN"))
dev.off()



#pseudotime_heatmap analysis
#######################################################
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

print(plot_pseudotime_heatmap(HSMM_myo[marker_gene,],
                              num_clusters = 5,
                              cores = 8,
                              show_rownames = T))


#TF_genes pseudotime_heatmap analysis
#######################################################
TF_genes <- c("NME2","EVI1","COMP1","CEBPB","SYNCRIP","GABP","TAXCREB","NPM1","CTR9","GTF2E2","ARHGAP35","YY1","TAF9B","HJURP","NEUROG3","SOX11","MAPK3","PSMB5","ZNF318","GTF2A2","FOXR2","PER1","ZZZ3")


TF_genes  <- intersect(rownames(HSMM_myo),TF_genes )


print(plot_pseudotime_heatmap(HSMM_myo[TF_genes,],
                              num_clusters = 3,
                              cores = 8,
                              show_rownames = T))

#plot_gene <- "IGHG4"


pdf(file=paste0("Figure2B_",Sys.Date(),"B_cell_plot_cell_trajectory_dbscan_puro1_",plot_gene,".pdf"))
plot_cell_trajectory(HSMM_myo, color_by =(matrix_temp_scale[plot_gene,]))+scale_colour_gradientn(colours = c("#2166ac","#67a9cf","#d1e5f0","#fddbc7","#ef8a62","#b2182b"),limits=c(min(matrix_temp_scale[plot_gene,]),max(matrix_temp_scale[plot_gene,])))+labs(title=paste(Cell_group,plot_gene))
dev.off()


pdf(file=paste0("Figure2B_",Sys.Date(),"B_cell_FeaturePlot_puro1_",plot_gene,".pdf"))
print(FeaturePlot(pbmc_part_metastatic_R, features = plot_gene,cols = c("#2166ac","#67a9cf","#d1e5f0","#fddbc7","#ef8a62","#b2182b"), pt.size=0.2))
dev.off()



#plot_gene <- "IGHG4"



pdf(file=paste0("Figure2B_",Sys.Date(),"B_cell_FeaturePlot_puro1_",plot_gene,".pdf"))
print(FeaturePlot(pbmc_pure_1, features = plot_gene,cols = c("#2166ac","#67a9cf","#d1e5f0","#fddbc7","#ef8a62","#b2182b") ))
dev.off()


print(FeaturePlot(pbmc, features = plot_gene,cols = c("#2166ac","#67a9cf","#d1e5f0","#fddbc7","#ef8a62","#b2182b") ))

pdf(file=paste0("Figure2B_",Sys.Date(),"B_cell_plot_cell_trajectory_dbscan_puro1_",plot_gene,".pdf"))
plot_cell_trajectory(HSMM_myo, color_by =(matrix_temp_scale[plot_gene,]))+scale_colour_gradientn(colours = c("#2166ac","#67a9cf","#d1e5f0","#fddbc7","#ef8a62","#b2182b"),limits=c(min(matrix_temp_scale[plot_gene,]),max(matrix_temp_scale[plot_gene,])))+labs(title=paste(Cell_group,plot_gene))
dev.off()






