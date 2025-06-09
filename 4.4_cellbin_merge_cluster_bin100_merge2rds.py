# coding: utf-8
import os,sys,argparse
del sys.path[4]
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.cluster import hierarchy
from scipy import sparse
from copy import deepcopy
sys.path.insert(0, '/jdfsbjcas1/ST_BJ/P21H28400N0232/guojingwen/software/others_code/gongchanghao')
from utils import getDefaultColors, removeBiasGenes
sys.path.insert(0, '/jdfsbjcas1/ST_BJ/P21H28400N0232/guojingwen/software/codes/BRCA_NCC_20220729')
from pytools import *

arg = argparse.ArgumentParser()
arg.add_argument('-a', '--adata_path', help='adata path') #/jdfsbjcas1/ST_BJ/P21H28400N0232/guojingwen/project/BRCA_NCC_20220729/ST/148TR_E1_16-1.2/Results/4_cellbin/7_new_cellbin/magic/all_hvgs/cellbin_annotated.h5ad
arg.add_argument('-m', '--mask_path', help='cell mask img', default='./') #/jdfsbjcas1/ST_BJ/P21H28400N0232/guojingwen/project/BRCA_NCC_20220729/ST/148BR_D1_13-1/Results/5_cellbin_230105/tissue_cut/SS200000148BR_D1_regist_mask_ft.tif
arg.add_argument('-o', '--out_path', help='output path', default='./') #/jdfsbjcas1/ST_BJ/P21H28400N0232/guojingwen/project/BRCA_NCC_20220729/ST/148BR_D1_13-1/Results/5_cellbin_230105/8.23_cellbin_merge_cluster_bin100_merge2rds
arg.add_argument('-b', '--batch', help='batch number')
args = arg.parse_args()

od = args.out_path
os.system(f"mkdir -p {od}")
os.chdir(od)
mask = args.mask_path
adata_path = args.adata_path
batch = args.batch

# 读取数据集
#adata = sc.read_h5ad("/jdfsbjcas1/ST_BJ/P21H28400N0232/guojingwen/project/BRCA_NCC_20220729/ST/148BR_D1_13-1.2/Results/4_cellbin/7_new_cellbin/SS200000148BR_D1.h5ad")
#anno = pd.read_table('/jdfsbjcas1/ST_BJ/P21H28400N0232/guojingwen/project/BRCA_NCC_20220729/ST/148BR_D1_13-1.2/Results/4_cellbin/7_new_cellbin/magic/change_hvgs/leiden/cellbin_metadata_annotated.txt', sep='\t')
#adata.obs = anno
#adata.write_h5ad('SS200000148BR_D1_raw_annotated_1.h5ad')
#adata = sc.read_h5ad("SS200000148BR_D1_raw_annotated_1.h5ad")
adata = sc.read_h5ad(adata_path) #注：h5ad转rds时，adata.X要为np.ndarray，不能为sparse matrix，因为matrix过大时sparse matrix会报错。在跑此脚本np.mean时，adata.X要为sparse matrix，因为np.ndarray跑出来不是三维的而是二维的。
adata = adata[adata.obs['spatial_domain'].isin([6,1,10,5,0])].copy()
adata = adata[adata.obs['batch'].isin([batch])].copy()
print(adata.obs['batch'].value_counts())
cluster_annotation = {6:'TS1',1:'TS2',10:'TS3',5:'TS4',0:'TS5',7:'TS6',2:'2',3:'3',4:'4',8:'8',9:'9',11:'11',12:'12',13:'13',14:'14'}
adata.obs['spatial_domain'] = adata.obs['spatial_domain'].map(cluster_annotation).astype('category')

os.system('ln -s %s %s'%(adata_path, od))
adata.obs.to_csv(od+"/"+adata_path.split('/')[-1].split('.')[0]+"_metadata.txt",sep = '\t')
if str(type(adata.X))=="<class 'scipy.sparse._csr.csr_matrix'>":
    adata.X = adata.X.A
adata_path2 = od+'/'+adata_path.rsplit('/',1)[1]
#adata.write_h5ad(adata_path2)
rds_path = adata_path2.rsplit('.',1)[0]+'.rds'
cmd_h5ad2rds = '/jdfsbjcas1/ST_BJ/P21H28400N0232/guojingwen/software/miniconda3/bin/Rscript /jdfsbjcas1/ST_BJ/P21H28400N0232/guojingwen/software/codes/BRCA_NCC_20220729/4.annh5ad2rds2.R --infile %s --outfile %s'%(adata_path2,rds_path)
print(cmd_h5ad2rds)
#os.system(cmd_h5ad2rds)
#os.rename(od+'/metadata.txt', od+'/'+adata_path2.split('/')[-1].split('.')[0]+'_metadata.txt')
#if str(type(adata.X))=="<class 'numpy.ndarray'>":
#    adata.X = sparse.csr_matrix(adata.X)


# 计算每个cluster的距离矩阵，根据距离阈值将聚类划分为簇，再合并簇的基因表达量
os.chdir(od)
adata_list = []
for category in adata.obs['spatial_domain'].unique():
    sub_adata = adata[adata.obs['spatial_domain'] == category].copy()
    adata_list.append(sub_adata)
merged_adatas = []
merged_cluster_labels = pd.DataFrame(columns=['cell_id', 'merged_cluster'])
adata_idx = 0
cell_num = 30
for adata in adata_list:
    print(adata_idx)
    # 计算距离矩阵
    dist_matrix = pdist(adata.obs[['x', 'y']])
    # 计算聚类树
    linkage = hierarchy.linkage(dist_matrix, method='ward')
    # 根据距离阈值将聚类划分为簇
    cluster_labels = hierarchy.fcluster(linkage, t=1000, criterion='distance')
    # 记录cell_id与新簇的对应关系
    merged_cluster_labels = pd.concat([merged_cluster_labels,
        pd.DataFrame(list(zip(adata.obs_names, ['%s_%s'%(i-1,adata_idx) for i in cluster_labels])), columns=['cell_id', 'merged_cluster'])])
    adata_idx += 1
    # 计算每个簇的中心坐标和gene count
    cluster_centers = []
    cluster_gene_counts = []
    if str(type(adata.X))=="<class 'numpy.ndarray'>":
        adata.X = sparse.csr_matrix(adata.X)
    for cluster_id in np.unique(cluster_labels):
        cluster_spots = adata.obs.index[cluster_labels == cluster_id]
        cluster_center = np.mean(adata.obs.loc[cluster_spots, ['x', 'y']], axis=0)
        cluster_gene_count = np.sum(adata[cluster_spots].X, axis=0)
        cluster_centers.append(cluster_center)
        cluster_gene_counts.append(cluster_gene_count)
    # 构建合并后的adata对象
    merged_adata = sc.AnnData(
        X=np.array(cluster_gene_counts)[:,0,:],
        obsm={'spatial': np.array(cluster_centers)},
        obs=pd.DataFrame(index=np.arange(len(cluster_centers))),
        var=adata.var,
    )
    merged_adata.obs['spatial_domain'] = adata.obs['spatial_domain'].iloc[0]
    merged_adatas.append(merged_adata)
merged_cluster_labels.to_csv(f'merged_cluster_labels_{batch}.txt', sep='\t', index=None)

# 合并所有合并后的adata对象
merged_adata = sc.AnnData.concatenate(*merged_adatas, join='outer', index_unique='_')
merged_adata

merged_adata.obs

# 保存结果
# cluster_number = merged_adata.obs['spatial_domain'].unique().shape[0]
# colors = getDefaultColors(cluster_number, type = 2)
# sc.pl.spatial(merged_adata, img_key="hires", color='spatial_domain',spot_size = 80,palette=colors,show=False)
# plt.savefig("merged_adata_spatial_cluster_refined_group.pdf",bbox_inches="tight")
merged_adata.var['mt'] = merged_adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(merged_adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
# sc.pl.violin(merged_adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
#              jitter=0.4, multi_panel=True,stripplot=False)
# plt.savefig("qc_violin.pdf",bbox_inches="tight")
merged_adata.obs['x'] = merged_adata.obsm['spatial'][:,0]
merged_adata.obs['y'] = merged_adata.obsm['spatial'][:,1]
sc.pp.normalize_total(merged_adata, inplace=True,target_sum=1e4)
sc.pp.log1p(merged_adata)
sc.pp.highly_variable_genes(merged_adata, flavor="seurat", n_top_genes=2000)
merged_adata = removeBiasGenes(merged_adata)
merged_adata.write_h5ad(f'adata_merge_{batch}.h5ad')

# sc.tl.rank_genes_groups(merged_adata, 'spatial_domain', method='wilcoxon', use_raw=False)
# plt.rcParams["figure.figsize"] = (7, 7)
# sc.pl.rank_genes_groups(merged_adata, n_genes=25, sharey=False)
# plt.savefig("rank_genes_groups.pdf")
# if str(type(merged_adata.X))=="<class 'numpy.ndarray'>":
#     merged_adata.X = sparse.csr_matrix(merged_adata.X)
# merged_adata.write_h5ad('adata_merge.h5ad')

# adata转成rds
merged_adata = sc.read_h5ad("adata_merge.h5ad")
os.system(f"export LD_LIBRARY_PATH=/jdfsbjcas1/ST_BJ/P21H28400N0232/guojingwen/software/miniconda3/lib:$LD_LIBRARY_PATH; /jdfsbjcas1/ST_BJ/P21H28400N0232/guojingwen/software/miniconda3/bin/Rscript /jdfsbjcas1/ST_BJ/P21H28400N0232/guojingwen/software/codes/BRCA_NCC_20220729/4.annh5ad2rds2.R --infile adata_merge.h5ad --outfile adata_merge.rds")



