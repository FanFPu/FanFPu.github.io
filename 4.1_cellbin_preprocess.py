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



data_path = sys.argv[1] ## gem.gz file
regist_path = sys.argv[2] ## regist tif
wd = os.path.dirname(os.path.abspath(data_path))
os.chdir(wd)
img_path = regist_path
python_path = '/jdfsbjcas1/ST_BJ/P21H28400N0232/gongchanghao/software/miniconda/envs/maskenv/bin/python'
cell_seg_api_path = '/jdfsbjcas1/ST_BJ/P21H28400N0232/guojingwen/project/BRCA_NCC_20220729/ST/148BR_D1_13-1.2/Results/4_cellbin/7_new_cellbin/cell_segmentation_v03_221209/cell_seg_api.pyc'
mask_path = os.path.basename(regist_path).replace("regist.tif","regist_mask.tif")
out_path = './cell_seg/' + mask_path
os.system("mkdir -p ./cell_seg/")
os.system("%s %s -i %s -o %s"%(python_path, cell_seg_api_path, img_path, out_path))


# tissue extraction
os.system("mkdir -p ./tissue_cut/")
from stereo.image.tissue_cut import SingleStrandDNATissueCut, DEEP
# Initial the TissueCut object
ssDNA_tissue_cut = SingleStrandDNATissueCut(
    seg_method=DEEP,
    src_img_path=regist_path,
    dst_img_path='./tissue_cut/',
    model_path="/jdfsbjcas1/ST_BJ/P21H28400N0232/guojingwen/project/BRCA_NCC_20220729/ST/148BR_D1_13-1.2/Results/4_cellbin/7_new_cellbin/weight_tissue_cut_tool_220304.hdf5"  # don't forget to download it, when using deep-learning method
)
# Real do the image transforming
ssDNA_tissue_cut.tissue_seg()
tissue_path = glob.glob("./tissue_cut/*tissue_cut.tif")[0]

# tissue extraction result filter
tissue_filtered_path = tissue_path.replace("tissue_cut.tif", "mask_ft.tif")
print("/jdfsbjcas1/ST_BJ/P21H28400N0232/gongchanghao/software/miniconda/envs/maskenv/bin/python /jdfsbjcas1/ST_BJ/P21H28400N0232/guojingwen/project/BRCA_NCC_20220729/ST/148BR_D1_13-1.2/Results/4_cellbin/7_new_cellbin/cell_segmentation_v03_221209/tissuecut_cellmask_filter.pyc -m %s -t %s -o ./tissue_cut"%(out_path, tissue_path))
os.system("/jdfsbjcas1/ST_BJ/P21H28400N0232/gongchanghao/software/miniconda/envs/maskenv/bin/python /jdfsbjcas1/ST_BJ/P21H28400N0232/guojingwen/project/BRCA_NCC_20220729/ST/148BR_D1_13-1.2/Results/4_cellbin/7_new_cellbin/cell_segmentation_v03_221209/tissuecut_cellmask_filter.pyc -m %s -t %s -o ./tissue_cut"%(out_path, tissue_path))

# cell correct
os.system("mkdir -p ./cell_correct_result")
os.system("/jdfsbjcas1/ST_BJ/P21H28400N0232/gongchanghao/software/miniconda/envs/maskenv/bin/python /jdfsbjcas1/ST_BJ/P21H28400N0232/guojingwen/project/BRCA_NCC_20220729/ST/148BR_D1_13-1.2/Results/4_cellbin/7_new_cellbin/cell_segmentation_v03_221209/GMMCorrect.pyc -m %s -g %s -o cell_correct_result -p 40 -t 20"%(tissue_filtered_path, data_path))
# os.system("/jdfsbjcas1/ST_BJ/P21H28400N0232/gongchanghao/software/miniconda/envs/maskenv/bin/python /jdfsbjcas1/ST_BJ/P21H28400N0232/gongchanghao/software/cell_segmentation_v03_221209/fast_little.pyc -m %s -g %s -o cell_correct_result"%(tissue_filtered_path, data_path))
# get result scanpy h5ad
corrected_gem = glob.glob("./cell_correct_result/data_adjust.txt")

data = st.io.read_gem(
        file_path= corrected_gem[0],
        sep='\t', bin_type="cell_bins",
        is_sparse=True)
out_h5ad = data_path.replace('gem.gz','h5ad')
out_rds = out_h5ad.replace('h5ad', 'rds')
adata = st.io.stereo_to_anndata(data,flavor='seurat',output=out_h5ad)
#os.system("export LD_LIBRARY_PATH=/jdfsbjcas1/ST_BJ/P21H28400N0232/gongchanghao/software/miniconda/envs/st/lib:$LD_LIBRARY_PATH; Rscript /jdfsbjcas1/ST_BJ/P21H28400N0232/gongchanghao/script/spatial_pipeline/bin/annh5ad2rds2.R --infile %s --outfile %s"%(out_h5ad, out_rds))





