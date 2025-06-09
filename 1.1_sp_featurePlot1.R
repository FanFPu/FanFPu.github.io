### Seurat clustering
library(Seurat)
library(dplyr)
library(cowplot)
library(patchwork)
#library(DoubletFinder)
library(ggplot2)
library(future)
library(Cairo)
options(bitmapType='cairo')
#library(psych)
#library(qgraph)
library(igraph)
library(Matrix)
#library(SeuratWrappers)
#library(pryr)
library(MuDataSeurat)

source("/jdfsbjcas1/ST_BJ/P21H28400N0232/wangdi/tumors_public/codes/scRNA_primary.R")
source("/jdfsbjcas1/ST_BJ/P21H28400N0232/wangdi/tumors_public/codes/colors.R")
source("/jdfsbjcas1/ST_BJ/P21H28400N0232/wangdi/tumors_public/codes/signatureGenes.R")
source("/jdfsbjcas1/ST_BJ/P21H28400N0232/wangdi/tumors_public/codes/0_initNEWanalysis.R")
source("/jdfsbjcas1/ST_BJ/P21H28400N0232/wangdi/tumors_public/codes/DEG_wilcox.R")

args <- commandArgs(T)
# topdir <- "/jdfsbjcas1/ST_BJ/P21H28400N0232/xieguixiang/BLCA/HBCP18/binbase/figures/"
workdir <- '/jdfsbjcas1/ST_BJ/P21H28400N0232/xieguixiang/BLCA/HBCP18/binbase/figures/score_plot2/'

if (!dir.exists(workdir)) {
  dir.create(workdir)
}

load('/jdfsbjcas1/ST_BJ/P21H28400N0232/xieguixiang/BLCA/HBCP18/binbase/seurat_spatialObj.RData')
setwd(workdir)

Basal <- list(c("KRT5","KRT14","KRT6A","CD44","KRT6","KRT17","CDH3","MMP14"))
Luminal <- list(c('KRT20','PPARG','FOXA1','GATA3','SNX31','UPK1A','UPK2','FGFR3','ESR1','PGR','KRT8','KRT18'))
EMT <- list(c("ABI3BP", "ACTA2", "ADAM12", "ANPEP", "APLP1", "AREG", "BASP1", "BDNF", "BGN", "BMP1", "CADM1", "CALD1", "CALU", "CAP2", "CAPG", "CCN1", "CCN2", "CD44", "CD59", "CDH11", "CDH2", "CDH6", "COL11A1", "COL12A1", "COL16A1", "COL1A1", "COL1A2", "COL3A1", "COL4A1", "COL4A2", "COL5A1", "COL5A2", "COL5A3", "COL6A2", "COL6A3", "COL7A1", "COL8A2", "COLGALT1", "COMP", "COPA", "CRLF1", "CTHRC1", "CXCL1", "CXCL12", "CXCL6", "CXCL8", "DAB2", "DCN", "DKK1", "DPYSL3", "DST", "ECM1", "ECM2", "EDIL3", "EFEMP2", "ELN", "EMP3", "ENO2", "FAP", "FAS", "FBLN1", "FBLN2", "FBLN5", "FBN1", "FBN2", "FERMT2", "FGF2", "FLNA", "FMOD", "FN1", "FOXC2", "FSTL1", "FSTL3", "FUCA1", "FZD8", "GADD45A", "GADD45B", "GAS1", "GEM", "GJA1", "GLIPR1", "GPC1", "GPX7", "GREM1", "HTRA1", "ID2", "IGFBP2", "IGFBP3", "IGFBP4", "IL15", "IL32", "IL6", "INHBA", "ITGA2", "ITGA5", "ITGAV", "ITGB1", "ITGB3", "ITGB5", "JUN", "LAMA1", "LAMA2", "LAMA3", "LAMC1", "LAMC2", "LGALS1", "LOX", "LOXL1", "LOXL2", "LRP1", "LRRC15", "LUM", "MAGEE1", "MATN2", "MATN3", "MCM7", "MEST", "MFAP5", "MGP", "MMP1", "MMP14", "MMP2", "MMP3", "MSX1", "MXRA5", "MYL9", "MYLK", "NID2", "NNMT", "NOTCH2", "NT5E", "NTM", "OXTR", "P3H1", "PCOLCE", "PCOLCE2", "PDGFRB", "PDLIM4", "PFN2", "PLAUR", "PLOD1", "PLOD2", "PLOD3", "PMEPA1", "PMP22", "POSTN", "PPIB", "PRRX1", "PRSS2", "PTHLH", "PTX3", "PVR", "QSOX1", "RGS4", "RHOB", "SAT1", "SCG2", "SDC1", "SDC4", "SERPINE1", "SERPINE2", "SERPINH1", "SFRP1", "SFRP4", "SGCB", "SGCD", "SGCG", "SLC6A8", "SLIT2", "SLIT3", "SNAI2", "SNTB1", "SPARC", "SPOCK1", "SPP1", "TAGLN", "TFPI2", "TGFB1", "TGFBI", "TGFBR3", "TGM2", "THBS1", "THBS2", "THY1", "TIMP1", "TIMP3", "TNC", "TNFAIP3", "TNFRSF11B", "TNFRSF12A", "TPM1", "TPM2", "TPM4", "VCAM1", "VCAN", "VEGFA", "VEGFC", "VIM", "WIPF1", "WNT5A"))



seurat_spatialObj <- AddModuleScore(seurat_spatialObj,
                          features = Basal,
                          ctrl = 100,
                          name = "Basal_score")
seurat_spatialObj <- AddModuleScore(seurat_spatialObj,
                          features = Luminal,
                          ctrl = 100,
                          name = "Luminal_score")
seurat_spatialObj <- AddModuleScore(seurat_spatialObj,
                          features = EMT,
                          ctrl = 100,
                          name = "EMT_score")

score_names <- c("Basal_score1", "Luminal_score1", "EMT_score1")

# 循环处理每个得分列
for (score_name in score_names) {
  # 计算当前得分列的 10% 最大值阈值
  threshold <- quantile(seurat_spatialObj@meta.data[[score_name]], 0.9)
  # 突出最大值的 10%
  seurat_spatialObj@meta.data[[score_name]] <- ifelse(seurat_spatialObj@meta.data[[score_name]] >= threshold, seurat_spatialObj@meta.data[[score_name]] * 2 , seurat_spatialObj@meta.data[[score_name]])
}
head(seurat_spatialObj@meta.data)

# 重命名列名
colnames(seurat_spatialObj@meta.data)[which(colnames(seurat_spatialObj@meta.data) == "Basal_score1")] <- "Basal"
colnames(seurat_spatialObj@meta.data)[which(colnames(seurat_spatialObj@meta.data) == "Luminal_score1")] <- "Luminal"
colnames(seurat_spatialObj@meta.data)[which(colnames(seurat_spatialObj@meta.data) == "EMT_score1")] <- "EMT"


pdf(file = paste0(workdir,"Basal_score.pdf"))
SpatialFeaturePlot(seurat_spatialObj, features = 'Basal', crop = TRUE, stroke=0)
dev.off()
pdf(file = paste0(workdir,"Luminal_score.pdf"))
SpatialFeaturePlot(seurat_spatialObj, features = 'Luminal', crop = TRUE, stroke=0)
dev.off()
pdf(file = paste0(workdir,"EMT_score.pdf"))
SpatialFeaturePlot(seurat_spatialObj, features = 'EMT', crop = TRUE, stroke=0)
dev.off()
