# -*- coding: utf-8 -*-
'''
Description: 
Author: gongch
Date: 2022-10-27 13:54:02
LastEditTime: 2022-10-27 13:54:03
LastEditors: gongchanghao
E-mail: gongchanghao@genomics.cn
'''

import os,sys
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import anndata
import random
import pandas as pd
from scvi.external import GIMVI
import scvi

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
            colors = colors[:n]
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

def removeBiasGenes(adata):
    malat1 = adata.var_names.str.startswith('MALAT1')
    MTgenes = adata.var_names.str.startswith('MT-')
    hb_genes = adata.var_names.str.contains('^HB[^(P)]')
    RPgenes = adata.var_names.str.startswith('RP') & adata.var_names.str.contains('-')
    RPgenes2 = adata.var_names.str.contains('^RP[SL]')
    CTCgenes = adata.var_names.str.startswith('CTC') & adata.var_names.str.contains('-')
    MIRgenes = adata.var_names.str.startswith('MIR')
    ACgenes = adata.var_names.str.contains('^AC[0-9]') & adata.var_names.str.contains('.')
    CTgenes = adata.var_names.str.startswith('CT') & adata.var_names.str.contains('-')
    LINCgenes = adata.var_names.str.contains('^LINC[0-9]')
    ALgenes = adata.var_names.str.contains('^AL') & adata.var_names.str.contains('.')
    remove_genes = malat1 | MTgenes | hb_genes | RPgenes | RPgenes2 | CTCgenes | MIRgenes | ACgenes | CTgenes | LINCgenes | ALgenes
    keep = np.invert(remove_genes)
    res = adata[:,keep]
    return res

args = sys.argv
indata = args[1]
mask = args[2]
outdir = args[3]
os.system("mkdir -p %s"%(outdir))
os.chdir(outdir)
st_adata = sc.read_h5ad(indata)
scvi.settings.num_threads = 50
# st_adata = sc.read_h5ad("/jdfsbjcas1/ST_BJ/P21H28400N0232/gongchanghao/project/PDAC/SC_reanalysis/result/10_cellbin_analysis/P2/cell_correct_result/SS200000148TR_D1.tissue_extraction.adjusted.cellbin.h5ad")
sc_adata = sc.read_h5ad("/jdfsbjcas1/ST_BJ/P21H28400N0232/guojingwen/project/BRCA_NCC_20220729/ST/41-1/Results/4_cellbin/4.2_scvi1/seurat_sc_P41.h5ad")
sc_adata.__dict__['_raw'].__dict__['_var'] = sc_adata.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})

#only use genes in both datasets
G = 3000

sc_adata.layers["counts"] = sc_adata.raw.X.astype(int).copy()
sc.pp.filter_cells(sc_adata, min_counts = 10)

st_adata.layers["counts"] = st_adata.X.astype(int).copy()
sc.pp.filter_cells(st_adata, min_counts = 10)


# filter genes to be the same on the spatial data
intersect = np.intersect1d(sc_adata.var_names, st_adata.var_names)
st_adata = st_adata[:, intersect].copy()
sc_adata = sc_adata[:, intersect].copy()
G = len(intersect)
print(G)

GIMVI.setup_anndata(st_adata,layer="counts")
GIMVI.setup_anndata(sc_adata, layer="counts", labels_key="cellTypes_new",batch_key="Patients")
model = GIMVI(sc_adata, st_adata)
model.train(max_epochs=200)
GIMVI.save(model, dir_path=outdir, overwrite=True, save_anndata=True)
model = GIMVI.load(outdir)
st_adata = sc.read_h5ad(outdir+'/adata_spatial.h5ad')
sc_adata = sc.read_h5ad(outdir+'/adata_seq.h5ad')
# st_adata2 = sc.read_h5ad("/jdfsbjcas1/ST_BJ/P21H28400N0232/gongchanghao/project/PDAC/spatial/cellbin/SS200000148TR_D1/cell_correct_result/SS200000148TR_D1.tissue_extraction.adjusted.cellbin.h5ad")
_, imputed = model.get_imputed_values(normalized=True)
st_adata.X = imputed
sc.pp.log1p(st_adata)
sc.pp.highly_variable_genes(st_adata, flavor="seurat", n_top_genes=2000)
sc.pp.pca(st_adata)
sc.pp.neighbors(st_adata)
sc.tl.umap(st_adata)
sc.tl.leiden(st_adata, key_added="res.0.4",resolution=0.4)
sc.tl.leiden(st_adata, key_added="res.0.8",resolution=0.8)
sc.tl.leiden(st_adata, key_added="res.1",resolution=1)
sc.tl.leiden(st_adata, key_added="res.1.4",resolution=1.4)
st_adata.write(outdir+"/cellbin_imputated.h5ad")


plt.rcParams["figure.figsize"] = (8, 8)
fig, axs = plt.subplots(2, 2, figsize=(16, 16))
sc.pl.umap(st_adata, color="res.0.4", wspace=0.4,show=False,ax=axs[0][0])
sc.pl.umap(st_adata, color="res.0.8", wspace=0.4,show=False,ax=axs[0][1])
sc.pl.umap(st_adata, color="res.1", wspace=0.4,show=False,ax=axs[1][0])
sc.pl.umap(st_adata, color="res.1.4", wspace=0.4,show=False,ax=axs[1][1])
plt.savefig(outdir+"/cluster_umap.pdf",bbox_inches="tight")
fig, axs = plt.subplots(2, 2, figsize=(25, 16))
cluster_number = st_adata.obs['res.0.4'].cat.categories.shape[0]
colors = getDefaultColors(cluster_number, type = 2)
sc.pl.spatial(st_adata, img_key="hires", color="res.0.4", size=1.5,spot_size = 20,palette=colors,show=False,ax=axs[0][0])
cluster_number = st_adata.obs['res.0.8'].cat.categories.shape[0]
colors = getDefaultColors(cluster_number, type = 2)
sc.pl.spatial(st_adata, img_key="hires", color="res.0.8", size=1.5,spot_size = 20,palette=colors,show=False,ax=axs[0][1])
cluster_number = st_adata.obs['res.1'].cat.categories.shape[0]
colors = getDefaultColors(cluster_number, type = 2)
sc.pl.spatial(st_adata, img_key="hires", color="res.1", size=1.5,spot_size = 20,palette=colors,show=False,ax=axs[1][0])
cluster_number = st_adata.obs['res.1.4'].cat.categories.shape[0]
colors = getDefaultColors(cluster_number, type = 2)
sc.pl.spatial(st_adata, img_key="hires", color="res.1.4", size=1.5,spot_size = 20,palette=colors,show=False,ax=axs[1][1])
plt.savefig("spatial_cluster_umap.pdf",bbox_inches="tight")


resolution = 'res.1' #choose best resolution
res = pd.DataFrame(st_adata.obs, columns = ["x", "y",resolution], index = st_adata.obs.index)
res.to_csv("bin1clu.txt",sep = '\t',index =False)
clusters = st_adata.obs[resolution].cat.categories
cluster_number = clusters.shape[0]
colors = getDefaultColors(cluster_number, type = 2)
flout = open("color.list",'w')
for i in range(cluster_number):
    if clusters[i] == 'low_quality':
        flout.write(clusters[i] + '\t#ffffff\n')
    else:
        flout.write(clusters[i] + '\t' + colors[i] + '\n')
flout.close()
os.system('/jdfsbjcas1/ST_BJ/P21H28400N0232/gongchanghao/bin/cell_bin_plot bin1clu.txt %s color.list cluster_plot.tif'%(mask))
os.system('mkdir -p cluster_split')
os.chdir("cluster_split")
flout = open("color.list",'w')
flout.write('1\t#ff0000\n')
flout.write('low_quality\t#ffffff\n')
flout.close()
for i in clusters:
    if i == 'low_quality':
        continue
    tmp = res.copy()
    tmp.loc[tmp[resolution] != i, [resolution]] = 'low_quality'
    tmp.loc[tmp[resolution] == i, [resolution]] = '1'
    tmp.to_csv("bin1clu_%s.txt"%(i),sep = '\t',index =False)
    os.system('/jdfsbjcas1/ST_BJ/P21H28400N0232/gongchanghao/bin/cell_bin_plot bin1clu_%s.txt %s color.list cluster_plot_%s.tif'%(i, mask, i))
st_adata = removeBiasGenes(st_adata)
if st_adata.raw:
    del st_adata.raw
if 'log1p' in st_adata.uns.keys():
    st_adata.uns['log1p']["base"] = None
sc.pp.highly_variable_genes(st_adata, flavor="seurat", n_top_genes=2000)
sc.tl.rank_genes_groups(st_adata, resolution, method='wilcoxon')
plt.rcParams["figure.figsize"] = (7,7)
sc.pl.rank_genes_groups(st_adata, n_genes=25, sharey=False)
plt.savefig("rank_genes_groups.pdf")

markers = pd.DataFrame(st_adata.uns['rank_genes_groups']['names']).head(3).stack().values.tolist()
markers = list(set(markers))
sc.pl.stacked_violin(st_adata, markers, groupby=resolution, rotation=90)
plt.savefig("marker_genes_violin.pdf")

degs = pd.DataFrame(st_adata.uns['rank_genes_groups']['names']).head(50)
degs.to_csv("degs.txt",sep = '\t')

'''
# get the latent representations for the sequencing and spatial data
latent_seq, latent_spatial = model.get_latent_representation()

#concatenate to one latent representation
latent_representation = np.concatenate([latent_seq, latent_spatial])
latent_adata = anndata.AnnData(latent_representation)

#labels which cells were from the sequencing dataset and which were from the spatial dataset
latent_labels = (['seq'] * latent_seq.shape[0]) + (['spatial'] * latent_spatial.shape[0])
latent_adata.obs['labels'] = latent_labels

#compute umap
sc.pp.neighbors(latent_adata, use_rep = 'X')
sc.tl.umap(latent_adata)
sc.tl.leiden(latent_adata, key_added="res.1",resolution=1)
sc.pl.umap(latent_adata, color="res.1", wspace=0.4)
plt.savefig("cluster_umap_latent.pdf")
sc.pl.umap(latent_adata, color = 'labels', show = True)
plt.savefig("labels_umap_latent.pdf")

cluster = latent_adata.obs['res.1'][sc_adata.shape[0]:]
cluster.index = st_adata.obs.index
st_adata.obs['scvi_cluster'] = cluster
clusters = st_adata.obs['scvi_cluster'].cat.categories
cluster_number = clusters.shape[0]
colors = getDefaultColors(cluster_number, type = 2)
sc.pl.spatial(st_adata, img_key="hires", color="scvi_cluster", size=1.5,spot_size = 20,palette=colors)
plt.savefig("spatial_cluster_umap_latent.pdf") 

# #save umap representations to original seq and spatial_datasets
sc_adata.obsm['X_umap'] = latent_adata.obsm['X_umap'][:sc_adata.shape[0]]
st_adata.obsm['X_umap'] = latent_adata.obsm['X_umap'][sc_adata.shape[0]:]

st_adata.write('scvi_cluster.h5ad')

# plot cellbin cluster
tmp = pd.merge(st_adata2.obs,st_adata.obs['scvi_cluster'],how = 'left', left_index = True, right_index = True)
tmp['scvi_cluster'] = tmp['scvi_cluster'].cat.add_categories("LQ").fillna('LQ')
st_adata2.obs = tmp
res = pd.DataFrame(st_adata2.obs, columns = ["x", "y","scvi_cluster"], index = st_adata2.obs.index)
res.to_csv("bin1clu_scvi.txt",sep = '\t',index =False)
clusters = st_adata2.obs['scvi_cluster'].cat.categories
cluster_number = clusters.shape[0]

colors = getDefaultColors(cluster_number, type = 1)
flout = open("color_scvi.list",'w')
for i in range(cluster_number):
    if clusters[i] == 'low_quality':
        flout.write(clusters[i] + '\t#ffffff\n')
    else:
        flout.write(clusters[i] + '\t' + colors[i] + '\n')
flout.close()
mask = '/jdfsbjcas1/ST_BJ/P21H28400N0232/gongchanghao/project/PDAC/SC_reanalysis/result/10_cellbin_analysis/P2/stereopy/SS200000148TR_D1_regist_mask.tif'
os.system('/jdfsbjcas1/ST_BJ/P21H28400N0232/gongchanghao/bin/cell_bin_plot bin1clu_merged.txt %s color.list cluster_plot_merged.tif'%(mask))
'''



