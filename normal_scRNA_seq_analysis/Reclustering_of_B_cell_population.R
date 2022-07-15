
library(Seurat)
library("monocle")
library(dplyr)
library(tidyr)
library(reshape2)



#################################################
###############set the working path##############
#Put ############################################
#/outs/filtered_feature_bc_matrix################
#and#############################################
#BandatB_cells_total_samples.txt#################
#and GSEA files##################################
#in the folder###################################
#################################################

setwd("work_folder")

################################################################
#########Clustering of all cell population###########
################################################################
#load files
X21_P.pbmc.data <- Read10X(data.dir = paste0("filtered_feature_bc_matrix"))
B_cells_cells_total<-read.delim("BandatB_cells_total_samples.txt")
REST<-read.delim("gsea_report_for_REST_1596397417231.tsv",header = T)
test_group<-read.delim("gsea_report_for_1_1596397417231.tsv",header = T)

#set out put folder under work_folder/
setwd("outputs/")

pbmc <- CreateSeuratObject(counts = X21_P.pbmc.data,names.field=2,names.delim="-" ,project = "S_aggr12", min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
plotVlnPlot<-VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc, features = rownames(pbmc))
pbmc <- ScaleData(pbmc, vars.to.regress = c("nFeature_RNA", "nCount_RNA", "percent.mt"), features = rownames(pbmc))

#Remove batch effects#
pbmc_list <- SplitObject(pbmc, split.by = "ident")
object.anchors <- FindIntegrationAnchors(object.list = pbmc_list, dims = 1:30)
pbmc <- IntegrateData(anchorset = object.anchors, dims = 1:30)
pbmc <- ScaleData(object = pbmc, verbose = FALSE)

pbmc <- RunPCA(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- RunUMAP(pbmc, dims = 1:10)
pbmc <- RunTSNE(pbmc, dims = 1:5)
pbmc <- FindClusters(pbmc, resolution = 0.2)
tsne_result<-data.frame(pbmc@reductions[["pca"]]@cell.embeddings,pbmc@reductions[["umap"]]@cell.embeddings,pbmc@reductions[["tsne"]]@cell.embeddings)
save(tsne_result,file="single_cell_cacer_atlas_cancertype_Seurat_individualgenes_tsne_result")
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, num.cores = 64)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)


####plot gene expression heatmap, distinguish those clusters as B cells#####
pdf(paste0("Seurat_DoHeatmap_",Sys.Date(),".pdf"))

print(DoHeatmap(pbmc, features = top10$gene) + NoLegend())

dev.off()

#plot
all_cells_total<-data.frame(Barcode=rownames(tsne_result),(data.frame(x=rownames(tsne_result)) %>% separate(x, c("A", "B"))))
sample<-recode(all_cells_total$B, "1"="S1","2"="S2","3"="S3","4"="S10","5"="S5","6"="S6","7"="S7","8"="S8","9"="S9","10"="S11","11"="S12","12"="S13","13"="S14")
all_cells_total<-data.frame(all_cells_total,sample)
meta<-recode(all_cells_total$sample, "S1"="no-metastasis","S2"="no-metastasis","S3"="liver-metastasis","S10"="liver-metastasis","S5"="liver-metastasis","S6"="no-metastasis","S7"="liver-metastasis","S8"="no-metastasis","S9"="liver-metastasis","S11"="no-metastasis","S12"="no-metastasis","S13"="liver-metastasis","S14"="no-metastasis")
side<-recode(all_cells_total$sample, "S1"="metastatic-R","S2"="metastatic-R","S3"="metastatic-R","S10"="metastatic-L","S5"="metastatic-R","S6"="metastatic-R","S7"="metastatic-R","S8"="metastatic-R","S9"="metastatic-L","S11"="metastatic-L","S12"="metastatic-L","S13"="metastatic-L","S14"="metastatic-L")
all_cells_list_tsne<-data.frame(tsne_result,all_cells_total,meta,side)
plot(match(all_cells_list_tsne$Barcode,rownames(all_cells_list_tsne)))

pdf(paste0("Figure_",Sys.Date(),"all_cells_list_umap_meta_sample.pdf"))
plot(all_cells_list_tsne$UMAP_1,all_cells_list_tsne$UMAP_2,col=as.factor(all_cells_list_tsne$meta), cex = 0.1)
dev.off()

################################################################
#########Re-clustering of B cell population###########
################################################################


X21_P.pbmc.data_B_cells_cells<-X21_P.pbmc.data[,which(X21_P.pbmc.data@Dimnames[[2]] %in% B_cells_cells_total$Barcode)]
pbmc_part <- CreateSeuratObject(counts = X21_P.pbmc.data_B_cells_cells,names.field=2,names.delim="-" ,project = "S_aggr12", min.cells = 3, min.features = 200)

pbmc_part@meta.data[["sample"]] <-B_cells_cells_total$sample
pbmc_part@meta.data[["meta"]] <-B_cells_cells_total$meta
pbmc_part@meta.data[["side"]] <-B_cells_cells_total$side
pbmc_part@meta.data[["Cluster"]] <-B_cells_cells_total$Cluster
pbmc_part@meta.data[["Cell_type"]] <-B_cells_cells_total$Cell_type
pbmc_part@meta.data[["LR_MN"]]<-(paste(pbmc_part@meta.data[["side"]],pbmc_part@meta.data[["meta"]],sep = "_"))

gc()

pbmc_part <- NormalizeData(pbmc_part, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc_part <- FindVariableFeatures(pbmc_part, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc_part)
pbmc_part <- ScaleData(pbmc_part, features = rownames(pbmc_part),block.size = 100,min.cells.to.block = 1000)

pbmc_part <- RunPCA(pbmc_part)
pbmc_part <- FindNeighbors(pbmc_part, dims = 1:10)
pbmc_part <- RunUMAP(pbmc_part, dims = 1:10)
pbmc_part <- RunTSNE(pbmc_part, dims = 1:10)
pbmc_part <- FindClusters(pbmc_part, resolution = 0.1)

pdf(file=paste0("Figure_",Sys.Date(),"B_cellreplot_tsne_4_groups.pdf"),10,8)

print(DimPlot(pbmc_part, reduction = "umap", group.by = "meta",cols =c("#7BC27A","#336087"),pt.size=1 ))
print(DimPlot(pbmc_part, reduction = "umap", group.by = "side",cols =c("#7BC27A","#336087"),pt.size=1))
print(DimPlot(pbmc_part, reduction = "tsne", group.by = "meta",cols =c("#7BC27A","#336087"),pt.size=1 ))
print(DimPlot(pbmc_part, reduction = "tsne", group.by = "side",cols =c("#7BC27A","#336087"),pt.size=1))

dev.off()


################################################################
#########Pseudotime analysis metastatic-L#############
################################################################

patient_name<-"S_aggr12" 
Selected_population <- "metastatic-L"

gc()


pbmc_part_metastatic_L <- pbmc_part[,which(pbmc_part$side ==Selected_population)]

matrix_temp <- pbmc_part_metastatic_L@assays[["RNA"]]@counts
matrix_temp_scale <- pbmc_part_metastatic_L@assays[["RNA"]]@scale.data

gc()

HSMM <- newCellDataSet(as.matrix(matrix_temp))

HSMM@phenoData@data[["side"]] <- pbmc_part_metastatic_L@meta.data[["side"]]
HSMM@phenoData@data[["meta"]] <- pbmc_part_metastatic_L@meta.data[["meta"]]
HSMM@phenoData@data[["sample"]] <- pbmc_part_metastatic_L@meta.data[["sample"]]
HSMM@phenoData@data[["Cell_type"]] <- pbmc_part_metastatic_L@meta.data[["Cell_type"]]
HSMM@phenoData@data[["Cluster"]] <- pbmc_part_metastatic_L@meta.data[["Cluster"]]

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


save(HSMM_myo,file=paste0("Pseudotime_monocle_obj",patient_name,Selected_population))

print(plot_cell_trajectory(HSMM_myo))
print(plot_cell_trajectory(HSMM_myo, color_by="Pseudotime"))
print(plot_cell_trajectory(HSMM_myo, color_by= "side")) 
print(plot_cell_trajectory(HSMM_myo, color_by= "meta"))
print(plot_cell_trajectory(HSMM_myo, color_by= "sample"))

print(plot_complex_cell_trajectory(HSMM_myo))
print(plot_complex_cell_trajectory(HSMM_myo, color_by="Pseudotime"))
print(plot_complex_cell_trajectory(HSMM_myo, color_by= "side"))
print(plot_complex_cell_trajectory(HSMM_myo, color_by= "meta"))
print(plot_complex_cell_trajectory(HSMM_myo, color_by= "sample"))

dev.off()


################################################################
########Pseudotime analysis metastatic-R#############
################################################################

patient_name<-"S_aggr12" 
Selected_population <- "metastatic-R"

gc()


pbmc_part_metastatic_R <- pbmc_part[,which(pbmc_part$side ==Selected_population)]

matrix_temp <- pbmc_part_metastatic_R@assays[["RNA"]]@counts
matrix_temp_scale <- pbmc_part_metastatic_R@assays[["RNA"]]@scale.data

gc()

HSMM <- newCellDataSet(as.matrix(matrix_temp))

HSMM@phenoData@data[["side"]] <- pbmc_part_metastatic_R@meta.data[["side"]]
HSMM@phenoData@data[["meta"]] <- pbmc_part_metastatic_R@meta.data[["meta"]]
HSMM@phenoData@data[["sample"]] <- pbmc_part_metastatic_R@meta.data[["sample"]]
HSMM@phenoData@data[["Cell_type"]] <- pbmc_part_metastatic_R@meta.data[["Cell_type"]]
HSMM@phenoData@data[["Cluster"]] <- pbmc_part_metastatic_R@meta.data[["Cluster"]]

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


save(HSMM_myo,file=paste0("Pseudotime_monocle_obj",patient_name,Selected_population))

print(plot_cell_trajectory(HSMM_myo))
print(plot_cell_trajectory(HSMM_myo, color_by="Pseudotime"))
print(plot_cell_trajectory(HSMM_myo, color_by= "side")) 
print(plot_cell_trajectory(HSMM_myo, color_by= "meta"))
print(plot_cell_trajectory(HSMM_myo, color_by= "sample"))

print(plot_complex_cell_trajectory(HSMM_myo))
print(plot_complex_cell_trajectory(HSMM_myo, color_by="Pseudotime"))
print(plot_complex_cell_trajectory(HSMM_myo, color_by= "side"))
print(plot_complex_cell_trajectory(HSMM_myo, color_by= "meta"))
print(plot_complex_cell_trajectory(HSMM_myo, color_by= "sample"))

dev.off()

################################################################
####Re-clustering of Activated B cell##########
################################################################

pbmc_part_actived_B_Cell <- pbmc_part[,which(pbmc_part$Cell_type == "actived_B_Cell")]

###############################
pbmc_part_actived_B_Cell@meta.data[["LR_MN"]]<-(paste(pbmc_part_actived_B_Cell@meta.data[["side"]],pbmc_part_actived_B_Cell@meta.data[["meta"]],sep = "_"))
gc()

pbmc_part_actived_B_Cell <- NormalizeData(pbmc_part_actived_B_Cell, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc_part_actived_B_Cell <- FindVariableFeatures(pbmc_part_actived_B_Cell, selection.method = "vst", nfeatures = 2000)
pbmc_part_actived_B_Cell <- ScaleData( pbmc_part_actived_B_Cell, vars.to.regress = c("percent.mt", "nCount_RNA"), features = rownames( pbmc_part_actived_B_Cell))
pbmc_part_actived_B_Cell <- RunPCA(pbmc_part_actived_B_Cell)
pbmc_part_actived_B_Cell <- FindNeighbors(pbmc_part_actived_B_Cell, dims = 1:10)
pbmc_part_actived_B_Cell <- RunUMAP(pbmc_part_actived_B_Cell, dims = 1:5)
pbmc_part_actived_B_Cell <- RunTSNE(pbmc_part_actived_B_Cell, dims = 1:5)
pbmc_part_actived_B_Cell <- FindClusters(pbmc_part_actived_B_Cell, resolution = 0.1)


pdf(file=paste0("Figure",patient_name,LSEC_population,Sys.Date(),"MT_meta_20200731_B_cellreplot_tsne_active_Bcell_groups.pdf"),10,9)
print(DimPlot(pbmc_part_actived_B_Cell, reduction = "umap", group.by = "meta",cols =c("#7BC27A","#336087"),pt.size=1, repel = T ))
dev.off()
pdf(file=paste0("Figure",patient_name,LSEC_population,Sys.Date(),"MT_side_20200731_B_cellreplot_tsne_active_Bcell_groups.pdf"),10,9)
print(DimPlot(pbmc_part_actived_B_Cell, reduction = "umap", group.by = "side",cols =c("#7BC27A","#336087"),pt.size=1, repel = T))
dev.off()
pdf(file=paste0("Figure",patient_name,LSEC_population,Sys.Date(),"MT_T_MATA_20200731_B_cellreplot_tsne_active_Bcell_groups.pdf"),10,9)
print(DimPlot(pbmc_part_actived_B_Cell, reduction = "tsne", group.by = "meta",cols =c("#7BC27A","#336087"),pt.size=1, repel = T ))
dev.off()
pdf(file=paste0("Figure",patient_name,LSEC_population,Sys.Date(),"MT_T_side_20200731_B_cellreplot_tsne_active_Bcell_groups.pdf"),10,9)
print(DimPlot(pbmc_part_actived_B_Cell, reduction = "tsne", group.by = "side",cols =c("#7BC27A","#336087"),pt.size=1, repel = T))
dev.off()
pdf(file=paste0("Figure",patient_name,LSEC_population,Sys.Date(),"MT_tsne_LR_MN_20200731_B_cellreplot_tsne_active_Bcell_groups.pdf"),10,7)
print(DimPlot(pbmc_part_actived_B_Cell, reduction = "tsne", group.by = "LR_MN",pt.size=1, repel = T))
dev.off()
pdf(file=paste0("Figure",patient_name,LSEC_population,Sys.Date(),"MT_LR_MN_20200731_B_cellreplot_tsne_active_Bcell_groups.pdf"),10,7)
print(DimPlot(pbmc_part_actived_B_Cell, reduction = "umap", group.by = "LR_MN",pt.size=1, repel = T))
dev.off()
pdf(file=paste0("Figure",patient_name,LSEC_population,Sys.Date(),"MT_tsne_cluseter_20200731_B_cellreplot_tsne_active_Bcell_groups.pdf"),10,10)
print(DimPlot(pbmc_part_actived_B_Cell, reduction = "tsne",pt.size=1, repel = T))
dev.off()

pdf(file=paste0("Figure",patient_name,LSEC_population,Sys.Date(),"MT_cluseter_20200731_B_cellreplot_tsne_active_Bcell_groups.pdf"),10,10)
print(DimPlot(pbmc_part_actived_B_Cell, reduction = "umap",pt.size=1, repel = T))
dev.off()


pdf(paste0("Figure_FeaturePlot",Sys.Date(),".pdf"),15,10)

print(FeaturePlot(pbmc_part_actived_B_Cell, reduction = "umap",features = c("IGHD","IGHM","TCL1A","LGALS1","IGHA1","IGHA2"),cols = c("grey","red"),ncol=3 ))

dev.off()

pbmc_part_actived_B_Cell.markers <- FindAllMarkers(pbmc_part_actived_B_Cell, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, num.cores = 64)
pbmc_part_actived_B_Cell.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
top10 <- pbmc_part_actived_B_Cell.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)


pdf(paste0("Figure_acvite_BC_DoHeatmap10_",Sys.Date(),".pdf"),4,5)

print(DoHeatmap(pbmc_part_actived_B_Cell, features = top10$gene) + NoLegend())

dev.off()


################################################################
##########Number of cells proportion bar plot#########
################################################################

seurat_clusters_LR_MN<-data.frame(seurat_clusters=pbmc_part_actived_B_Cell$seurat_clusters,LR_MN=pbmc_part_actived_B_Cell$LR_MN)
seurat_clusters_LR_MN<-data.frame(t(table(seurat_clusters_LR_MN)))
barpt<-ggplot(melt(seurat_clusters_LR_MN), aes(fill=seurat_clusters, y=value, x=LR_MN)) + 
  geom_bar(position="fill", stat="identity")+
  coord_flip()

pdf(paste0("Figure_bar_plot_",Sys.Date(),".pdf"),15,10)
barpt
dev.off()

############################################################################
#############make GSEA files and use GSEA software to do GSEA analysis######
############################################################################

##############################################################
######################Figure 2E plot GSEA#####################
##############################################################
#load GSEA output 


REST_20<-subset(REST, GS.DETAILS!="")
REST_20$Group<-c("REST")
test_group_20<-subset(test_group, GS.DETAILS!="")
test_group_20$Group<-c("Cluster_1")
All<-rbind(REST_20, test_group_20)

bar_all_GSEA<-ggplot(All, aes(reorder(NAME, NES), NES, fill=Group)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values=c("#00BFC4", "black"))+
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Pathways NES from GSEA (Activated B cells)")

pdf(file="1596397417231Pathways_NES_from_GSEA.pdf",width=10, height=10)
bar_all_GSEA
dev.off()

