library("Seurat")
#library(biomaRt)
library(dplyr)
library(tidyr)
library("ggplot2")
library(randomcoloR)
library(RColorBrewer)
require("GSVA")

setwd("E:/XXXXXXXX/S_aggr12/S_aggr12/14sample_fil")

gc()
# all cells analysis
pbmc_allcell<-readRDS("20210228_14sample_fil.rds")

# add information
all_cells_total<-data.frame(Barcode=colnames(pbmc_allcell),(data.frame(x=colnames(pbmc_allcell)) %>% separate(x, c("A", "B"))))
sample<-recode(all_cells_total$B, "1"="S1","2"="S2","3"="S3","4"="S10","5"="S5","6"="S6","7"="S7","8"="S8","9"="S9","10"="S11","11"="S12","12"="S13","13"="S14")
pbmc_allcell@meta.data[["sample"]] <-sample
cell_type<-recode(pbmc_allcell$seurat_clusters , "7" = "Fibroblasts", "12" = "Endothelial_cells","0"="CD4_T_cells","1"="Cytotoxic_T_cells","2"="Activated/Plasma_like_B_cell","3"="KRT18_tumor_cells","6"="KRT8_tumor_cells","8"="Proliferative_tumor_cells","4"="Activated_like_B_cell","5"="Macrophage","9"="Mast_cells","10"="Plasma_cells","11"="Plasma_like_B_cells")
pbmc_allcell@meta.data[["Cell_type"]]<-  cell_type


single_cell_sample_means <-  data.frame(matrix(nrow=25455,ncol=0))

for (sample_num in unique(pbmc_allcell$sample)) 
  {
  
  
  #sample_num = "S1"
  
  pbmc_tem <- subset(pbmc_allcell, sample == sample_num )
  
  pbmc_tem <- as.matrix(pbmc_tem@assays[["RNA"]]@data) 
  
  
  single_cell_sample_means[,sample_num]<- (rowMeans(pbmc_tem))
  
  
}


rownames(single_cell_sample_means) <- names(rowMeans(pbmc_tem))


Heavy_type_list <-c("IGHG4","IGHG1","IGHG3","IGHV","IGHG2","JSRP1","DNAAF1")#Heavy_type_list

Light_type_list<- c("IGLC3","IGLC2","IGLV1","IGLV6","IGLC7","DERL3")#Light_type_list





Heavy_type_GSVA<-try(t(gsva(as.matrix(single_cell_sample_means),list(Light_type_list))))

Light_type_GSVA<-try(t(gsva(as.matrix(single_cell_sample_means),list(Heavy_type_list))))

HeavytoLight_Score <- Heavy_type_GSVA/Light_type_GSVA #Save Organize and Add Comment Columns

HeavytoLight_Score_gsva<-read.delim(file="HeavytoLight_Score_gsva.txt")

res.wilcox.test <- wilcox.test(HeavytoLight_Score ~  meta_or_not , data = HeavytoLight_Score_gsva,alternative = c("less"))


# library
library(ggplot2)
library(car)
# grouped boxplot
HeavytoLight_Score_gsvaplot<-ggplot(HeavytoLight_Score_gsva, aes(x=meta_or_not, y=HeavytoLight_Score, fill=side)) + 
  geom_boxplot()+
  scale_fill_brewer(palette="Dark2")


pdf(file="HeavytoLight_Score_gsva_all_single_cell.pdf",width=10, height=10)
print(HeavytoLight_Score_gsvaplot)
dev.off()