import os, sys
del sys.path[4]
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
import warnings

warnings.filterwarnings("ignore")

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["font.serif"] = ["Arial"]
sc.settings.set_figure_params(
    dpi=100, dpi_save=300, frameon=False, facecolor="white", fontsize=16, vector_friendly=True, figsize=(5, 5)
)
sc._settings.ScanpyConfig(figdir="./", n_jobs=30)

od = "/jdfsbjcas1/ST_BJ/P21H28400N0232/gongchanghao/project/PDAC/spatial/cellbin_v3/2.9_cellbin_merge_metacell"
os.system(f"mkdir -p {od}")
os.chdir(od)

from sklearn.cluster import BisectingKMeans

st_adata = sc.read(
    "/jdfsbjcas1/ST_BJ/P21H28400N0232/gongchanghao/project/PDAC/spatial/cellbin_v3/1.1_cellbin_rawdata/filtered_data/merged_adata5.h5ad"
)

# 根据分类标签合成大细胞
def merge_big_cell_v2(adata, resolution, prefix, n=30):
    adata_list = []
    for category in adata.obs[resolution].unique():
        sub_adata = adata[adata.obs[resolution] == category].copy()
        adata_list.append(sub_adata)
    merged_adatas = []
    merged_cluster_labels = pd.DataFrame(columns=["cell_id", "merged_cluster"])
    for i, st_adata in enumerate(adata_list):
        t = st_adata.shape[0] // n
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
        adata_cluster.obs["celltype"] = st_adata.obs["celltype"].iloc[0]
        adata_cluster.obs["cellsubtype"] = st_adata.obs["cellsubtype"].iloc[0]
        adata_cluster.obs["region"] = st_adata.obs["region"].iloc[0]
        merged_adatas.append(adata_cluster)
    merged_cluster_labels.to_csv(f"merged_cluster_labels_{prefix}.txt", sep="\t", index=None)
    # 合并所有合并后的adata对象
    merged_adata = sc.AnnData.concatenate(*merged_adatas, join="outer")
    return merged_adata


from scipy.sparse import csr_matrix
from STutils.pl import getDefaultColors

merged_adata_list = []
for i, sample in enumerate(samplelist):
    adata = st_adata[st_adata.obs["batch"] == sample, :].copy()
    count_dict = adata.obs[merge_label].value_counts().to_dict()
    less_than_100 = [k for k, v in count_dict.items() if v < 100]
    less_index = adata.obs[merge_label].isin(less_than_100)
    adata = adata[~less_index]
    adata.obs[merge_label] = adata.obs[merge_label].cat.remove_unused_categories()
    merged_adata = merge_big_cell_v2(adata, merge_label, samplelist[i], 20)
    sparse_X = csr_matrix(merged_adata.X)
    merged_adata.X = sparse_X
    merged_adata.obs["x"] = merged_adata.obsm["spatial"][:, 0]
    merged_adata.obs["y"] = merged_adata.obsm["spatial"][:, 1]
    # merged_adata = merged_adata_list[i]
    merged_adata.write_h5ad(f"merged_adata_{sample}.h5ad")
    merged_adata_list.append(merged_adata)
    merged_cluster_labels = pd.read_csv(f"merged_cluster_labels_{sample}.txt", sep="\t", index_col=0)
    adata.obs["merged_cluster_labels"] = merged_cluster_labels
    adata.obs["merged_cluster_labels"][adata.obs["merged_cluster_labels"].isna()] = "Unassigned_Cell"
    adata.obs["merged_cluster_labels"] = adata.obs["merged_cluster_labels"].astype("category")
    mask = f"/jdfsbjcas1/ST_BJ/P21H28400N0232/gongchanghao/project/PDAC/spatial/rawdata/{sample}/{sample}_regist_mask_ft.tif"
    res = pd.DataFrame(adata.obs, columns=["x", "y", "merged_cluster_labels"], index=adata.obs.index)
    res.to_csv(f"bin1clu_{sample_dict[sample]}.txt", sep="\t", index=False)

    clusters = adata.obs["merged_cluster_labels"].cat.categories
    cluster_number = clusters.shape[0]
    colors = getDefaultColors(cluster_number, type=2)
    flout = open(f"color.list_{sample_dict[sample]}", "w")
    for i in range(cluster_number):
        flout.write(clusters[i] + "\t" + colors[i] + "\n")
    flout.close()
    os.system(
        f"/jdfsbjcas1/ST_BJ/P21H28400N0232/gongchanghao/bin/cell_bin_plot bin1clu_{sample_dict[sample]}.txt {mask} color.list_{sample_dict[sample]} cluster_plot_{sample_dict[sample]}.tif"
    )

all_merged_adata = sc.AnnData.concatenate(*merged_adata_list, join="outer", batch_categories=samplelist)
all_merged_adata.obs["sample"] = all_merged_adata.obs["batch"].map(sample_dict)
all_merged_adata.write("all_merged_adata.h5ad")