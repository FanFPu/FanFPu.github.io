library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
options(stringsAsFactors = FALSE, future.globals.maxSize=20000*1024^2)

od <- '/jdfsbjcas1/ST_BJ/P21H28400N0232/wangdi/xieguixiang/BRCAFigures/SPACEL/region_cellchat/merge_IBCAR1_IBCAR2/P13P41/'
setwd(od)

x <- load('/jdfsbjcas1/ST_BJ/P21H28400N0232/wangdi/xieguixiang/BRCAFigures/SPACEL/region_cellchat/IBCAR1/P13P41/cellchat_interaction.RData')
cellchat_high_TIL <- get(x)
y <- load('/jdfsbjcas1/ST_BJ/P21H28400N0232/wangdi/xieguixiang/BRCAFigures/SPACEL/region_cellchat/IBCAR2/P13P41/cellchat_interaction.RData')
cellchat_low_TIL <- get(y)

cellchat <- mergeCellChat(list(cellchat_high_TIL, cellchat_low_TIL), add.names = c("IBCAR1", "IBCAR2"), cell.prefix = TRUE)
pdf(paste0(od, "compare_number.pdf"),width=3,height=3)
compareInteractions(cellchat, show.legend = F, group = c(1,2))
compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
dev.off()

pdf(paste0(od, "compare_number_bycell.pdf"),width=6,height=4)
netVisual_heatmap(cellchat)
netVisual_heatmap(cellchat, measure = "weight")
dev.off()

pdf(paste0(od, "rankNET.pdf"),width=4,height=7)
rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
dev.off()

library(ComplexHeatmap)
cellchat_high_TIL <- netAnalysis_computeCentrality(cellchat_high_TIL, slot.name = "netP")
cellchat_low_TIL <- netAnalysis_computeCentrality(cellchat_low_TIL, slot.name = "netP")
object.list <- list(IBCAR1 = cellchat_high_TIL, IBCAR2 = cellchat_low_TIL)
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 16)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 16)
pdf(paste0(od, "compare_outgoing_signal_bycell.pdf"),width=8,height=13)
draw(ht1 + ht2, ht_gap = unit(1.5, "cm"))
dev.off()
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 16, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 16, color.heatmap = "GnBu")
pdf(paste0(od, "compare_incoming_signal_bycell.pdf"),width=8,height=13)
draw(ht1 + ht2, ht_gap = unit(1.5, "cm"))
dev.off()
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 16, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 16, color.heatmap = "OrRd")
pdf(paste0(od, "compare_all_signal_bycell.pdf"),width=8,height=13)
draw(ht1 + ht2, ht_gap = unit(1.5, "cm"))
dev.off()

##### Compare L-R pairs
#pdf(paste0(od, "L-R_InvEpi_Immune_compare.pdf"),width=6,height=12)
#netVisual_bubble(cellchat, sources.use = c('Epithelial_Invasive_cancer','myCAF'), targets.use = c('B_cell','Macrophage','Plasmocyte','T_cell'),  comparison = c(1, 2), angle.x = 45)
#dev.off()

pdf(paste0(od, "IncreaseorDecrease_IBC_fib_signal_IBCAR2.pdf"),width=15,height=11)
gg1 <- netVisual_bubble(cellchat, sources.use = c('IBC cells','Fibroblasts'), targets.use = c('IBC cells','Fibroblasts'),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in IBCAR2 samples", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = c('IBC cells','Fibroblasts'), targets.use = c('IBC cells','Fibroblasts'),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in IBCAR2 samples", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2
dev.off()

# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "IBCAR2"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = "IBCAR2",ligand.logFC = 0.2, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = "IBCAR1",ligand.logFC = -0.1, receptor.logFC = -0.1)
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)
pairLR.use.up = net.up[, "interaction_name", drop = F]
print(pairLR.use.up)
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = c('Fibroblasts'), targets.use = c('IBC cells','Fibroblasts','Macrophages','T cells','Pericytes','Endothelial cells'),
                        comparison = c(1, 2),  angle.x = 45, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
#gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = c('IBC cells','Fibroblasts'), targets.use = c('Fibroblasts'), 
                        #comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
pdf(paste0(od, "Up-regulated_fib_sender_signal.pdf"),width=8,height=35)
gg1
dev.off()


#gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = c('IBC cells'), targets.use = c('Macrophages','Endothelial cells','Plasmocytes','Pericytes'),
                        #comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = c('IBC cells'), targets.use = c('IBC cells','Macrophages','Endothelial cells','T cells','B cells'), 
                        comparison = c(1, 2),  angle.x = 45, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
pdf(paste0(od, "Down-regulated_IBC_sender_signal.pdf"),width=5,height=12)
gg2
dev.off()

gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = c('T cells'), targets.use = c('IBC cells'), 
                        comparison = c(1, 2),  angle.x = 45, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
pdf(paste0(od, "Down-regulated_Tcell_sender_signal.pdf"),width=3,height=8)
gg2
dev.off()

gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = c('IBC cells'), targets.use = c('IBC cells','Macrophages','Endothelial cells','T cells','B cells','Fibroblasts'), 
                        comparison = c(1, 2),  angle.x = 45, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
pdf(paste0(od, "Up-regulated_IBC_sender_signal.pdf"),width=6,height=28)
gg2
dev.off()

gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = c('Pericytes'), targets.use = c('IBC cells'), 
                        comparison = c(1, 2),  angle.x = 45, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
pdf(paste0(od, "Down-regulated_pericytes_sender_signal.pdf"),width=4,height=13)
gg2
dev.off()


##### plot gene expression
cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("IBCAR1", "IBCAR2")) # set factor level
pdf(paste0(od, "GeneExpression.pdf"),width=13,height=100)
plotGeneExpression(cellchat, signaling = "Other", split.by = "datasets", colors.ggplot = T)
dev.off()