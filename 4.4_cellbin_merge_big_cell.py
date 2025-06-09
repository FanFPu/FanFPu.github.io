
# %%
import os, sys
# del sys.path[4]
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import BisectingKMeans
from scipy.sparse import csr_matrix
sys.path.insert(0, '/jdfsbjcas1/ST_BJ/P21H28400N0232/guojingwen/software/others_code/gongchanghao')
#from utils import getDefaultColors, removeBiasGenes
def removeBiasGenes(adata):
    IGgenes = adata.var_names.str.startswith('IG')
    COLgenes = adata.var_names.str.startswith('COL')
    malat1 = adata.var_names.str.startswith('MALAT1')
    MTgenes = adata.var_names.str.startswith('MT')
    hb_genes = adata.var_names.str.contains('^HB[^(P)]')
    RPgenes = adata.var_names.str.startswith('RP') & adata.var_names.str.contains('-')
    RPgenes2 = adata.var_names.str.contains('^RP[SL]')
    CTCgenes = adata.var_names.str.startswith('CTC') & adata.var_names.str.contains('-')
    MIRgenes = adata.var_names.str.startswith('MIR')
    ACgenes = adata.var_names.str.contains('^AC[0-9]') & adata.var_names.str.contains('.')
    CTgenes = adata.var_names.str.startswith('CT') & adata.var_names.str.contains('-')
    LINCgenes = adata.var_names.str.contains('^LINC[0-9]')
    ALgenes = adata.var_names.str.contains('^AL') & adata.var_names.str.contains('.')

    MESgenes = adata.var_names.isin(['ACTA2','ACTG2','ACTG1','TAGLN','DES','ACTC1','LUM','BGN','DCN','POSTN','MRC1',"CALD1","TPM1","MGP","C1S","PDPN",
                                      'FLT1','PLVAP','SPARCL1','PECAM1','CD31','HSPG2','VWF','CDH5','SELE','VCAM1','ENG','PDGFRB'])
    Tcellgenes = adata.var_names.isin(['TRAC','SPN','TAGAP','IL7R','PTPRC','IL10RA','LTB','RGS1','TRBC2','LCP2','FCMR','IL2RB','NLRC3','IL21R','IL18R1','SLAMF1','LAT2','CD2','CD52','KLRB1',
                                        "CD3D", "CD3E", "CD3G","CD8A", "CD8B", "GZMK", "CD4","TNFRSF4"])
    Plasmocytegenes = adata.var_names.isin(['FAM30A','MZB1','FCRL5','JCHAIN','LAX1','THEMIS2','SLAMF7','LY9','JSRP1','NCKAP1L'])
    Macrogenes = adata.var_names.isin(['LILRB5','MPEG1','MS4A7','FCGR2A','LYVE1','SIGLEC1','SIRPB2','RGS1','THEMIS2','ITGAX','C1QA','KCNE1','TYROBP','CSF3R','C1QB','CNR2','ADGRE1',
                                        'CD163','SLC11A1','APOC1','FCER1G','FCGR3A'])
    Bcellgenes = adata.var_names.isin(['MS4A1','BLK','CD79B','TAGAP','TNFRSF13C','FCMR','P2RX5','TLR10','FCRL1','CYBB','SCIMP','IRF8','LILRB1','THEMIS2','CR2','FCRL2','FCRL5','CD19','CD21','CD79A','CD79B','BLNK'])
    MYLgenes = adata.var_names.str.startswith('MYL')
    remove_genes = IGgenes | COLgenes | malat1 | MTgenes | hb_genes | RPgenes | RPgenes2 | CTCgenes | MIRgenes | ACgenes | CTgenes | LINCgenes | ALgenes | MESgenes | Tcellgenes | Plasmocytegenes | Macrogenes | Bcellgenes | MYLgenes
    # remove_genes = malat1 | MTgenes | hb_genes | RPgenes | RPgenes2 | CTCgenes | MIRgenes | ACgenes | CTgenes | LINCgenes | ALgenes
    keep = np.invert(remove_genes)
    res = adata[:,keep]
    return res

# 根据分类标签合成大细胞
def merge_big_cell(adata, resolution, n=30):
    adata_list = []
    for category in adata.obs[resolution].unique():
        sub_adata = adata[adata.obs[resolution] == category].copy()
        adata_list.append(sub_adata)
    merged_adatas = []
    merged_cluster_labels = pd.DataFrame(columns=["cell_id", "merged_cluster"])
    for i, st_adata in enumerate(adata_list):
        t = st_adata.shape[0] // n
        if t == 0:
            t = 1
        X = st_adata.obs[["x", "y"]].values
        spectral_clustering = BisectingKMeans(n_clusters=t, bisecting_strategy="largest_cluster", random_state=0)
        cluster_labels = spectral_clustering.fit_predict(X)
        # cluster_sizes = np.bincount(cluster_labels)
        # sorted_clusters = np.argsort(cluster_sizes)[::-1]
        # max_cluster_index = sorted_clusters[0]
        # while max_cluster_index
        cluster = st_adata.obs[resolution].iloc[0]
        adata_idx = adata.obs[resolution].unique().to_list().index(cluster)
        merged_cluster_labels_i = pd.DataFrame(
            list(zip(st_adata.obs_names, ["%s_%s" % (i - 1, adata_idx) for i in cluster_labels])),
            columns=["cell_id", "merged_cluster"],
        )
        merged_cluster_labels = pd.concat([merged_cluster_labels, merged_cluster_labels_i], axis=0)
        # 计算每个簇的中心坐标和gene count
        cluster_centers = []
        cluster_gene_counts = []
        cluster_cell_counts = []
        for cluster_id in np.unique(cluster_labels):
            cluster_spots = st_adata.obs.index[cluster_labels == cluster_id]
            cluster_center = np.mean(st_adata.obs.loc[cluster_spots, ["x", "y"]], axis=0)
            cluster_gene_count = np.sum(st_adata[cluster_spots].X, axis=0)
            cluster_centers.append(cluster_center)
            cluster_gene_counts.append(cluster_gene_count)
            cluster_cell_counts.append(cluster_spots.shape[0])
        # 构建合并后的adata对象
        adata_cluster = sc.AnnData(
            X=np.array(cluster_gene_counts)[:, 0, :],
            obsm={"spatial": np.array(cluster_centers)},
            obs=pd.DataFrame(
                cluster_cell_counts, index=np.arange(len(cluster_centers)), columns=["cluster_cell_counts"]
            ),
            var=adata.var,
        )
        adata_cluster.obs["merged_cluster"] = st_adata.obs[resolution].iloc[0]
        adata_cluster.obs["spatial_domain"] = st_adata.obs["spatial_domain"].iloc[0]
        # adata_cluster.obs["cellsubtype"] = st_adata.obs["cellsubtype"].iloc[0]
        # adata_cluster.obs["region"] = st_adata.obs["region"].iloc[0]
        merged_adatas.append(adata_cluster)
    merged_cluster_labels.to_csv("merged_cluster_labels.txt", sep="\t", index=None)
    # 合并所有合并后的adata对象
    merged_adata = sc.AnnData.concatenate(*merged_adatas, join="outer")
    return merged_adata

# samplelist = [
#     "HBCP1", "HBCP3A", "HBCP3B", "HBCP3C","HBCP5A","HBCP5B","HBCP6","HBCP7","HBCP8","HBCP9","HBCP10","HBCP11","HBCP12","HBCP13","HBCP14","HBCP15","HBCP17","HBCP18"]
samplelist = [
    "HBCP3", "HBCP18"]
# samplelist = [
#     "HBCP3A", "HBCP3B", "HBCP3C","HBCP6","HBCP18"]
# for sample in samplelist:
od = "/jdfsbjcas1/ST_BJ/P21H28400N0232/xieguixiang/BLCA/ST/immune/with_fibsubtype/"
os.system("mkdir -p %s"%(od))
os.chdir(od)

for sample in samplelist:    
    adata = sc.read_h5ad(f'/jdfsbjcas1/ST_BJ/P21H28400N0232/xieguixiang/BLCA/ST/immune/with_fibsubtype/{sample}.h5ad')
    adata.obs['annotated_cluster'] = adata.obs['annotated_cluster'].astype('str')
    # adata = adata[adata.obs['annotated_cluster'].isin(['MIBC','Metastatic_tumor','NMIBC'])].copy()
    if 'annotated_cluster' in adata.obs.columns and 'immune_subtype' in adata.obs.columns:
        for idx in adata.obs.index:
            if adata.obs['annotated_cluster'][idx] in ['Bcell','DCS','Macrophage','Mastcell','Plasmocyte','Tcell']:
                adata.obs.at[idx, 'annotated_cluster'] = adata.obs['immune_subtype'][idx]
    if 'annotated_cluster' in adata.obs.columns and 'fibsubtype' in adata.obs.columns:
        for idx in adata.obs.index:
            if adata.obs['annotated_cluster'][idx] == 'Fibroblast':
                adata.obs.at[idx, 'annotated_cluster'] = adata.obs['fibsubtype'][idx]
    adata.obs['annotated_cluster'] = adata.obs['annotated_cluster'].astype('category')
    # adata = adata[adata.obs['spatial_domain'].isin([11,2,10,0,3,14])].copy()
    # adata = adata[adata.obs['batch'].isin(["HBCP18"])].copy()
    # adata = adata[adata.obs['Region'].isin(["Red","Green","Blue"])].copy()
    # adata.obs['Region'] = adata.obs['Region'].astype('str')
    # adata.obs.loc[adata.obs["Region"]=='Red', ["Region"]]="ET"
    # adata.obs.loc[adata.obs["Region"]=='Green', ["Region"]]="ST"
    # adata.obs.loc[adata.obs["Region"]=='Blue', ["Region"]]="IT"
    # adata.obs['Region'] = adata.obs['Region'].astype('category')
    # adata.obs.loc[adata.obs['batch'].isin(['HBCP1','HBCP7','HBCP9','HBCP13','HBCP14','HBCP15','HBCP17','HBCP5A','HBCP5B']), ["Region"]]="NMIBC"
    # adata.obs.loc[adata.obs['batch'].isin(['HBCP3B','HBCP3C','HBCP6','HBCP8','HBCP10','HBCP3A','HBCP12','HBCP18']), ["Region"]]="MIBC"
    # cluster_annotation = {6:'TS1',1:'TS2',10:'TS3',5:'TS4',0:'TS4',7:'7',2:'2',3:'3',4:'4',8:'8',9:'9',11:'11',12:'12',13:'13',14:'14'}
    cluster_annotation = {0:'SD0',1:'SD1',2:'SD2',3:'SD3',4:'SD4',5:'SD5',6:'SD6',7:'SD7',8:'SD8',9:'SD9',10:'SD10',11:'SD11',12:'SD12',13:'SD13',14:'SD14'}
    adata.obs['spatial_domain'] = adata.obs['spatial_domain'].map(cluster_annotation).astype('category')

    # print(adata.obs["spatial_domain"].value_counts())
    # print(adata.obs["Region"].value_counts())
    # print(adata.obs["batch"].value_counts())
    # adata.write_h5ad('Total_tumor_smallCell.h5ad')
    #adata = adata[adata.obs["annotated_cluster"].isin(['T_cell','B_cell','Macrophage','Plasmocyte','Epithelial_Invasive_cancer','myCAF'])].copy()
    #adata = adata[adata.obs["annotated_cluster"]!="low_quality"].copy()
    st_ad_list = []
    for sd in np.unique(adata.obs['spatial_domain'].values):
        adata_select = adata[adata.obs['spatial_domain'] == sd]
        merged_adata = merge_big_cell(adata_select, 'annotated_cluster', 20)
        sparse_X = csr_matrix(merged_adata.X)
        merged_adata.X = sparse_X
        merged_adata.obs['x'] = merged_adata.obsm['spatial'][:, 0]
        merged_adata.obs['y'] = merged_adata.obsm['spatial'][:, 1]
        st_ad_list.append(merged_adata)
    # merged_adata = removeBiasGenes(merged_adata)
    adata = sc.AnnData.concatenate(*st_ad_list, join='outer')
    
    adata.write_h5ad(f'{sample}_bigcell.h5ad')



