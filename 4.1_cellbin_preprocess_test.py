# -*- coding: utf-8 -*-
'''
Description: 
Author: gongch
Date: 2022-12-12 13:53:57
LastEditTime: 2022-12-12 13:53:58
LastEditors: gongchanghao
E-mail: gongchanghao@genomics.cn
'''

#%%
import warnings

warnings.filterwarnings('ignore')
import stereo as st
import os,sys
import glob
import h5py
import numpy as np
import pandas as pd
from anndata import AnnData
from typing_extensions import Literal
from scipy.sparse import csr_matrix
from shapely.geometry import Point, MultiPoint

from stereo.io import h5ad
from stereo.core.cell import Cell
from stereo.core.gene import Gene
from stereo.core.stereo_exp_data import StereoExpData
from stereo.utils.read_write_utils import ReadWriteUtils
from stereo.log_manager import logger

data_path = sys.argv[1] ## gem.gz file
regist_path = sys.argv[2] ## regist tif
wd = os.path.dirname(os.path.abspath(data_path))
os.chdir(wd)
img_path = regist_path
python_path = '/jdfsbjcas1/ST_BJ/P21H28400N0232/gongchanghao/software/miniconda/envs/maskenv/bin/python'
cell_seg_api_path = '/jdfsbjcas1/ST_BJ/P21H28400N0232/guojingwen/project/BRCA_NCC_20220729/ST/148BR_D1_13-1.2/Results/4_cellbin/7_new_cellbin/cell_segmentation_v03_221209/cell_seg_api.pyc'
mask_path = os.path.basename(regist_path).replace("regist.tif","regist_mask.tif")
out_path = './cell_seg/' + mask_path

# get result scanpy h5ad
corrected_gem = glob.glob("./cell_correct_result/data_adjust.txt")

data = st.io.read_gem(
        file_path= corrected_gem[0],
        sep='\t', bin_type="cell_bins",
        is_sparse=True)
out_h5ad = data_path.replace('gem.gz','h5ad')
out_rds = out_h5ad.replace('h5ad', 'rds')

def stereo_to_anndata(data=data, flavor='seurat', sample_id="sample", reindex=False, output=out_h5ad, split_batches=True
):
    """
    Transform the StereoExpData object into Anndata format.  

    Parameters
    -----------------------
    data
        the input StereoExpData object.
    flavor
        if you want to convert the output file into h5ad of Seurat, please set `'seurat'`.
    sample_id
        the sample name which will be set as `orig.ident` in obs.
    reindex
        if `True`, the cell index will be reindexed as `{sample_id}:{position_x}_{position_y}` format.
    output
        the path to output h5ad file.
	split_batches
		Whether to save each batch to a single file if it is a merged data, default to True.
    Returns
    -----------------
    An object of Anndata.
    """
    if data.merged and split_batches:
        from os import path
        from ..utils.data_helper import split
        data_list = split(data)
        batch = np.unique(data.cells.batch)
        adata_list = []
        if output is not None:
            name, ext = path.splitext(output)
        for bno, d in zip(batch, data_list):
            if output is not None:
                boutput = f"{name}-{d.sn}{ext}"
            else:
                boutput = None
            adata = stereo_to_anndata(d, flavor=flavor, sample_id=sample_id, reindex=reindex, output=output,
                                      split_batches=False)
            adata_list.append(adata)
        return adata_list

    from scipy.sparse import issparse

    #if data.tl.raw is None:
        #logger.error('convert to AnnData should have raw data')
        #raise Exception

    exp = data.exp_matrix if issparse(data.exp_matrix) else csr_matrix(data.exp_matrix)

    cells = data.cells.to_df()
    cells.dropna(axis=1, how='all', inplace=True)
    print(cells)
    genes = data.genes.to_df()
    genes.dropna(axis=1, how='all', inplace=True)

    adata = AnnData(X=exp, dtype=np.float64, obs=cells, var=genes)
    #adata.raw = AnnData(X=data.tl.raw.exp_matrix, dtype=np.float64, obs=data.tl.raw.cells.to_df(),
                        #var=data.tl.raw.genes.to_df())

    ##sample id
    logger.info(f"Adding {sample_id} in adata.obs['orig.ident'].")
    adata.obs['orig.ident'] = pd.Categorical([sample_id] * adata.obs.shape[0], categories=[sample_id])
    if data.position is not None:
        logger.info(f"Adding data.position as adata.obsm['spatial'] .")
        adata.obsm['spatial'] = data.position
        # adata.obsm['X_spatial'] = data.position
        logger.info(f"Adding data.position as adata.obs['x'] and adata.obs['y'] .")
        adata.obs['x'] = pd.DataFrame(data.position[:, 0], index=data.cell_names.astype('str'))
        adata.obs['y'] = pd.DataFrame(data.position[:, 1], index=data.cell_names.astype('str'))

    if data.sn is not None:
        if isinstance(data.sn, str):
            sn_list = [['-1', data.sn]]
        else:
            sn_list = []
            for bno, sn in data.sn.items():
                sn_list.append([bno, sn])
        adata.uns['sn'] = pd.DataFrame(sn_list, columns=['batch', 'sn'])

    for key in data.tl.key_record.keys():
        if len(data.tl.key_record[key]) > 0:
            if key == 'hvg':
                res_key = data.tl.key_record[key][-1]
                logger.info(f"Adding data.tl.result['{res_key}'] into adata.var .")
                adata.uns[key] = {'params': {}, 'source': 'stereopy', 'method': key}
                for i in data.tl.result[res_key]:
                    if i == 'mean_bin':
                        continue
                    adata.var[i] = data.tl.result[res_key][i]
            elif key == 'sct':
                res_key = data.tl.key_record[key][-1]
                # adata.uns[res_key] = {}
                logger.info(f"Adding data.tl.result['{res_key}'] into adata.uns['sct_'] .")
                adata.uns['sct_counts'] = csr_matrix(data.tl.result[res_key][0]['counts'].T)
                adata.uns['sct_data'] = csr_matrix(data.tl.result[res_key][0]['data'].T)
                adata.uns['sct_scale'] = csr_matrix(data.tl.result[res_key][0]['scale.data'].T.to_numpy())
                adata.uns['sct_scale_genename'] = list(data.tl.result[res_key][0]['scale.data'].index)
                adata.uns['sct_top_features'] = list(data.tl.result[res_key][1]['top_features'])
                adata.uns['sct_cellname'] = list(data.tl.result[res_key][1]['umi_cells'].astype('str'))
                adata.uns['sct_genename'] = list(data.tl.result[res_key][1]['umi_genes'])
            elif key in ['pca', 'umap', 'tsne']:
                # pca :we do not keep variance and PCs(for varm which will be into feature.finding in pca of seurat.)
                res_key = data.tl.key_record[key][-1]
                sc_key = f'X_{key}'
                logger.info(f"Adding data.tl.result['{res_key}'] into adata.obsm['{sc_key}'] .")
                adata.obsm[sc_key] = data.tl.result[res_key].values
            elif key == 'neighbors':
                # neighbor :seurat use uns for conversion to @graph slot, but scanpy canceled neighbors of uns at present.
                # so this part could not be converted into seurat straightly.
                for res_key in data.tl.key_record[key]:
                    sc_con = 'connectivities' if res_key == 'neighbors' else f'{res_key}_connectivities'
                    sc_dis = 'distances' if res_key == 'neighbors' else f'{res_key}_distances'
                    logger.info(f"Adding data.tl.result['{res_key}']['connectivities'] into adata.obsp['{sc_con}'] .")
                    logger.info(f"Adding data.tl.result['{res_key}']['nn_dist'] into adata.obsp['{sc_dis}'] .")
                    adata.obsp[sc_con] = data.tl.result[res_key]['connectivities']
                    adata.obsp[sc_dis] = data.tl.result[res_key]['nn_dist']
                    logger.info(f"Adding info into adata.uns['{res_key}'].")
                    adata.uns[res_key] = {}
                    adata.uns[res_key]['connectivities_key'] = sc_con
                    adata.uns[res_key]['distance_key'] = sc_dis
                    # adata.uns[res_key]['connectivities'] = data.tl.result[res_key]['connectivities']
                    # adata.uns[res_key]['distances'] = data.tl.result[res_key]['nn_dist']
            elif key == 'cluster':
                for res_key in data.tl.key_record[key]:
                    logger.info(f"Adding data.tl.result['{res_key}'] into adata.obs['{res_key}'] .")
                    adata.obs[res_key] = pd.DataFrame(data.tl.result[res_key]['group'].values,
                                                      index=data.cells.cell_name.astype('str'))
            elif key in ('gene_exp_cluster', 'cell_cell_communication'):
                for res_key in data.tl.key_record[key]:
                    logger.info(f"Adding data.tl.result['{res_key}'] into adata.uns['{key}@{res_key}']")
                    adata.uns[f"{key}@{res_key}"] = data.tl.result[res_key]
            elif key == 'regulatory_network_inference':
                for res_key in data.tl.key_record[key]:
                    logger.info(f"Adding data.tl.result['{res_key}'] in adata.uns['{res_key}'] .")
                    regulon_key = f'{res_key}_regulons'
                    adata.uns[regulon_key] = data.tl.result[res_key]['regulons']
                    auc_matrix_key = f'{res_key}_auc_matrix'
                    adata.uns[auc_matrix_key] = data.tl.result[res_key]['auc_matrix']
                    adjacencies_key = f'{res_key}_adjacencies'
                    adata.uns[adjacencies_key] = data.tl.result[res_key]['adjacencies']
            else:
                continue

    if data.tl.raw is not None:
        if flavor == 'seurat':
            # keep same shape between @counts and @data for seurat,because somtimes dim of sct are not the same.
            logger.info(f"Adding data.tl.raw.exp_matrix as adata.uns['raw_counts'] .")
            adata.uns['raw_counts'] = data.tl.raw.exp_matrix if issparse(data.tl.raw.exp_matrix) \
                else csr_matrix(data.tl.raw.exp_matrix)
            adata.uns['raw_cellname'] = list(data.tl.raw.cell_names.astype(str))
            adata.uns['raw_genename'] = list(data.tl.raw.gene_names)
            if data.tl.raw.position is not None and reindex:
                logger.info(f"Reindex as adata.uns['raw_cellname'] .")
                raw_sample = pd.DataFrame(['sample'] * data.tl.raw.cell_names.shape[0],
                                          index=data.tl.raw.cell_names.astype('str'))
                raw_x = pd.DataFrame(data.tl.raw.position[:, 0].astype(str), index=data.tl.raw.cell_names.astype('str'))
                raw_y = pd.DataFrame(data.tl.raw.position[:, 1].astype(str), index=data.tl.raw.cell_names.astype('str'))
                new_ix = np.array(raw_sample + "_" + raw_x + "_" + raw_y).tolist()
                adata.uns['raw_cellname'] = new_ix
        else:
            logger.info(f"Adding data.tl.raw.exp_matrix as adata.raw .")
            raw_exp = data.tl.raw.exp_matrix
            raw_genes = data.tl.raw.genes.to_df()
            raw_genes.dropna(axis=1, how='all', inplace=True)
            raw_adata = AnnData(X=raw_exp, var=raw_genes, dtype=np.float64)
            adata.raw = raw_adata

    if reindex:
        logger.info(f"Reindex adata.X .")
        new_ix = (adata.obs['orig.ident'].astype(str) + ":" + adata.obs['x'].astype(str) + "_" +
                  adata.obs['y'].astype(str)).tolist()
        adata.obs.index = new_ix
        if 'sct_cellname' in adata.uns.keys():
            logger.info(f"Reindex as adata.uns['sct_cellname'] .")
            adata.uns['sct_cellname'] = new_ix
    
    adata.obs['orig.ident'].astype(str)
    adata.obs['x'].astype(str)
    adata.obs['y'].astype(str)
    del adata.obs['cell_point']
    print(adata.obs)

    if flavor == 'seurat':
        logger.info(f"Rename QC info.")
        adata.obs.rename(columns={'total_counts': "nCount_Spatial", "n_genes_by_counts": "nFeature_Spatial",
                                  "pct_counts_mt": 'percent.mito'}, inplace=True)
        # if 'X_pca' not in list(adata.obsm.keys()):
        # logger.info(f"Creating fake info. Please ignore X_ignore in your data.")
        # adata.obsm['X_ignore'] = np.zeros((adata.obs.shape[0], 2))

    logger.info(f"Adding data.attr in adata.uns.")
    if data.offset_x is not None:
        adata.uns['offset_x'] = data.offset_x
    if data.offset_y is not None:
        adata.uns['offset_y'] = data.offset_y
    if data.attr is not None:
        for key, value in data.attr.items():
            adata.uns[key] = value

    logger.info(f"Finished conversion to anndata.")

    if output is not None:
        adata.write_h5ad(output)
        logger.info(f"Finished output to {output}")

    return adata


adata = stereo_to_anndata(data,flavor='seurat',output=out_h5ad)
#os.system("export LD_LIBRARY_PATH=/jdfsbjcas1/ST_BJ/P21H28400N0232/gongchanghao/software/miniconda/envs/st/lib:$LD_LIBRARY_PATH; Rscript /jdfsbjcas1/ST_BJ/P21H28400N0232/gongchanghao/script/spatial_pipeline/bin/annh5ad2rds2.R --infile %s --outfile %s"%(out_h5ad, out_rds))





