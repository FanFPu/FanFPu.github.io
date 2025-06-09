#加载分析包
current_lib_paths <- .libPaths()
new_lib_paths <- c(current_lib_paths[3], current_lib_paths[-3])
.libPaths(new_lib_paths)
.libPaths()

library(nichenetr)
library(tidyverse)
#library(circlize)
#library(Seurat)
#library(Matrix)
#library(Cairo)
#library(MuDataSeurat)
#options(bitmapType = "cairo")
library(DiagrammeR)
library(DiagrammeRsvg)
setwd('E:/BLCA/nichenet')


## 路径分析主要用到前两个文件；后3个文件主要用于数据来源/类型注释
weighted_networks = readRDS("weighted_networks_nsga2r_final.rds")
ligand_tf_matrix = readRDS("ligand_target_matrix_nsga2r_final.rds")

lr_network = readRDS("lr_network_human_allInfo_30112033.rds")
sig_network = readRDS("signaling_network_human_21122021.rds")
gr_network = readRDS("gr_network_human_21122021.rds")

## 交代感兴趣的配体与受体
ligands_all = c("HMGB1") # 可以是多个
targets_all = c("IFI16")
## 分析/提取路径
active_signaling_network = get_ligand_signaling_path(
  ligand_tf_matrix = ligand_tf_matrix, 
  ligands_all = ligands_all, 
  targets_all = targets_all, 
  weighted_networks = weighted_networks)
lapply(active_signaling_network, head)
## 权重标准化
active_signaling_network_min_max = active_signaling_network
active_signaling_network_min_max$sig = active_signaling_network_min_max$sig %>% 
  mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)
active_signaling_network_min_max$gr = active_signaling_network_min_max$gr %>% 
  mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)


## 转为dgr_graph对象
graph_min_max = diagrammer_format_signaling_graph(
  signaling_graph_list = active_signaling_network_min_max, 
  ligands_all = ligands_all, 
  targets_all = targets_all, 
  sig_color = "indianred", 
  gr_color = "steelblue")


## 可视化
DiagrammeR::render_graph(graph_min_max, 
                         output="graph", # graph, visNetwork
                         layout = "nicely" # nicely, circle, tree, kk, and fr
)
