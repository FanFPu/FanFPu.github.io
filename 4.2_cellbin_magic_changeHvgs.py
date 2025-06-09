# -*- coding: utf-8 -*-
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
import scanpy.external as sce
import os,sys
import random
#from pytools import *

def getDefaultColors(n, type = 1):
    if type == 1:
        colors = ["#ff1a1a", "#1aff1a", "#1a1aff", "#ffff1a", "#ff1aff",
                "#ff8d1a", "#7cd5c8", "#c49a3f", "#5d8d9c", "#90353b",
                "#507d41", "#502e71", "#1B9E77", "#c5383c", "#0081d1",
                "#674c2a", "#c8b693", "#aed688", "#f6a97a", "#c6a5cc",
                "#798234", "#6b42c8", "#cf4c8b", "#666666", "#ffd900",
                "#feb308", "#cb7c77", "#68d359", "#6a7dc9", "#c9d73d"]
    elif type == 2:
        if n <= 14:
            colors = ["#437BFE", "#FEC643", "#43FE69", "#FE6943", "#C643FE",
                  "#43D9FE", "#B87A3D", "#679966", "#993333", "#7F6699",
                  "#E78AC3", "#333399", "#A6D854", "#E5C494"]
        elif n <= 20:
            colors = ["#87b3d4", "#d5492f", "#6bd155", "#683ec2", "#c9d754",
                  "#d04dc7", "#81d8ae", "#d34a76", "#607d3a", "#6d76cb",
                  "#ce9d3f", "#81357a", "#d3c3a4", "#3c2f5a", "#b96f49",
                  "#4e857e", "#6e282c", "#d293c8", "#393a2a", "#997579"]
        elif n <= 30:
            colors = ["#628bac", "#ceda3f", "#7e39c9", "#72d852", "#d849cc",
                  "#5e8f37", "#5956c8", "#cfa53f", "#392766", "#c7da8b",
                  "#8d378c", "#68d9a3", "#dd3e34", "#8ed4d5", "#d84787",
                  "#498770", "#c581d3", "#d27333", "#6680cb", "#83662e",
                  "#cab7da", "#364627", "#d16263", "#2d384d", "#e0b495",
                  "#4b272a", "#919071", "#7b3860", "#843028", "#bb7d91"]
        else:
            colors = ["#982f29", "#5ddb53", "#8b35d6", "#a9e047", "#4836be",
                  "#e0dc33", "#d248d5", "#61a338", "#9765e5", "#69df96",
                  "#7f3095", "#d0d56a", "#371c6b", "#cfa738", "#5066d1",
                  "#e08930", "#6a8bd3", "#da4f1e", "#83e6d6", "#df4341",
                  "#6ebad4", "#e34c75", "#50975f", "#d548a4", "#badb97",
                  "#b377cf", "#899140", "#564d8b", "#ddb67f", "#292344",
                  "#d0cdb8", "#421b28", "#5eae99", "#a03259", "#406024",
                  "#e598d7", "#343b20", "#bbb5d9", "#975223", "#576e8b",
                  "#d97f5e", "#253e44", "#de959b", "#417265", "#712b5b",
                  "#8c6d30", "#a56c95", "#5f3121", "#8f846e", "#8f5b5c"]
    elif type == 3:
        colors = ["#588dd5", "#c05050", "#07a2a4", "#f5994e",
                "#9a7fd1", "#59678c", "#c9ab00", "#7eb00a"]
    elif type == 4:
        colors = ["#FC8D62", "#66C2A5", "#8DA0CB", "#E78AC3",
                "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3"]    
    elif type == 5:
        colors = ["#c14089", "#6f5553", "#E5C494", "#738f4c",
                "#bb6240", "#66C2A5", "#2dfd29", "#0c0fdc"]
    if n:
        if n <= len(colors):
            colors = colors[0:n]
        else:
            step = 16777200 // (n - len(colors)) - 2
            add_colors = []
            tmp = random.sample(range(step),1)[0]
            for i in range(n-len(colors)):
                hextmp = str(hex(tmp)).replace("0x",'')
                if len(hextmp) == 5:
                    hextmp = "0"+hextmp
                add_colors.append("#" + hextmp)
                tmp = tmp + step
            colors = colors + add_colors
    return colors

args = sys.argv
indata = args[1]
#scdata = args[2]
mask = args[2]
outdir = args[3]
os.system("mkdir -p %s"%(outdir))
os.system("mkdir -p %s/{cluster_umap,spatial_cluster_umap,rank_genes_groups,marker_genes_violin}"%(outdir))
os.chdir(outdir)
adata = sc.read_h5ad(indata)
adata.var_names_make_unique()
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
fig, axs = plt.subplots(1, 4, figsize=(15, 4))
sns.distplot(adata.obs["total_counts"], kde=False, ax=axs[0])
sns.distplot(adata.obs["total_counts"][adata.obs["total_counts"] < 10000], kde=False, bins=40, ax=axs[1])
sns.distplot(adata.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[2])
sns.distplot(adata.obs["n_genes_by_counts"][adata.obs["n_genes_by_counts"] < 4000], kde=False, bins=60, ax=axs[3])
plt.savefig("preQC.pdf")
# sc.pp.filter_cells(adata, min_counts=5000)
# sc.pp.filter_cells(adata, max_counts=35000)
# adata = adata[adata.obs["pct_counts_mt"] < 2]
sc.pp.filter_genes(adata, min_cells=10)
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
sce.pp.magic(adata, name_list='all_genes', knn=5)
adata.obs.to_csv("cellbin_metadata_magic.txt",sep = '\t')
adata.write("cellbin_magic.h5ad")

### change highly variable genes
sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)
#scdata = sc.read_h5ad(scdata)
#adata.var['highly_variable'] = False
#print(sum(adata.var.highly_variable))
#for i in scdata.var.index:
  #print(i)
  #if i in adata.var.index:
    #adata.var.loc[i, 'highly_variable'] = True
#print(sum(adata.var.highly_variable))

sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata, key_added="res.0.2",resolution=0.2)
sc.tl.leiden(adata, key_added="res.0.4",resolution=0.4)
sc.tl.leiden(adata, key_added="res.0.6",resolution=0.6)
sc.tl.leiden(adata, key_added="res.0.8",resolution=0.8)
sc.tl.leiden(adata, key_added="res.1",resolution=1)
sc.tl.leiden(adata, key_added="res.1.2",resolution=1.2)
adata.obs.to_csv("cellbin_metadata_Ownhvgs.txt",sep = '\t')
adata.write("cellbin_Ownhvgs.h5ad")

resolution = 'res.0.6'
res = pd.DataFrame(adata.obs, columns = ["x", "y",resolution], index = adata.obs.index)
res.to_csv("bin1clu_res06.txt",sep = '\t',index =False)
clusters = adata.obs[resolution].cat.categories
cluster_number = clusters.shape[0]
colors = getDefaultColors(cluster_number, type = 2)
flout = open("color_res06.list",'w')
for i in range(cluster_number):
    if clusters[i] == 'low_quality':
        flout.write(clusters[i] + '\t#ffffff\n')
    else:
        flout.write(clusters[i] + '\t' + colors[i] + '\n')
flout.close()
os.system('/jdfsbjcas1/ST_BJ/P21H28400N0232/gongchanghao/bin/cell_bin_plot bin1clu_res06.txt %s color_res06.list cluster_plot_res06.tif'%(mask))

resolution = 'res.0.8'
res = pd.DataFrame(adata.obs, columns = ["x", "y",resolution], index = adata.obs.index)
res.to_csv("bin1clu_res08.txt",sep = '\t',index =False)
clusters = adata.obs[resolution].cat.categories
cluster_number = clusters.shape[0]
colors = getDefaultColors(cluster_number, type = 2)
flout = open("color_res08.list",'w')
for i in range(cluster_number):
    if clusters[i] == 'low_quality':
        flout.write(clusters[i] + '\t#ffffff\n')
    else:
        flout.write(clusters[i] + '\t' + colors[i] + '\n')
flout.close()
os.system('/jdfsbjcas1/ST_BJ/P21H28400N0232/gongchanghao/bin/cell_bin_plot bin1clu_res08.txt %s color_res08.list cluster_plot_res08.tif'%(mask))

resolution = 'res.1'
res = pd.DataFrame(adata.obs, columns = ["x", "y",resolution], index = adata.obs.index)
res.to_csv("bin1clu_res1.txt",sep = '\t',index =False)
clusters = adata.obs[resolution].cat.categories
cluster_number = clusters.shape[0]
colors = getDefaultColors(cluster_number, type = 2)
flout = open("color_res1.list",'w')
for i in range(cluster_number):
    if clusters[i] == 'low_quality':
        flout.write(clusters[i] + '\t#ffffff\n')
    else:
        flout.write(clusters[i] + '\t' + colors[i] + '\n')
flout.close()
os.system('/jdfsbjcas1/ST_BJ/P21H28400N0232/gongchanghao/bin/cell_bin_plot bin1clu_res1.txt %s color_res1.list cluster_plot_res1.tif'%(mask))


