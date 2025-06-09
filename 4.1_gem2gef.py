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

#data = st.io.read_gef(file_path='/jdfsbjcas1/ST_BJ/P21H28400N0232/guojingwen/project/BRCA_NCC_20220729/ST/29-1/Results/4_cellbin/cellbin/cell_correct_result/SS200000149TL_E5.tissue_extraction.adjusted.cellbin.gef', bin_type='cell_bins')
#print(data.cells['area'])
data_path = sys.argv[1] ## gem.gz file
wd = os.path.dirname(os.path.abspath(data_path))
os.chdir(wd)
corrected_gem = glob.glob("./cell_correct_result/data_adjust.txt")
data = st.io.read_gem(
        file_path= corrected_gem[0],
        sep='\t', bin_type="cell_bins",
        is_sparse=True)

def write_mid_gef(data: StereoExpData, output: str):
    """
    Write the StereoExpData object into a GEF (.h5) file.

    Parameters
    ---------------------
    data
        the input StereoExpData object.
    output
        the path to output file.

    Returns
    ---------------------
    None
    """
    logger.info("The output standard gef file only contains one expression matrix with mid count."
                "Please make sure the expression matrix of StereoExpData object is mid count without normaliztion.")
    import numpy.lib.recfunctions as rfn
    final_exp = []  # [(x_1,y_1,umi_1),(x_2,y_2,umi_2)]
    final_gene = []  # [(A,offset,count)]
    exp_np = data.exp_matrix.toarray()

    for i in range(exp_np.shape[1]):
        gene_exp = exp_np[:, i]
        c_idx = np.nonzero(gene_exp)[0]  # idx for all cells
        zipped = np.concatenate((data.position[c_idx], gene_exp[c_idx].reshape(c_idx.shape[0], 1)), axis=1)
        for k in zipped:
            final_exp.append(k)

        # count
        g_len = len(final_gene)
        last_offset = 0 if g_len == 0 else final_gene[g_len - 1][1]
        last_count = 0 if g_len == 0 else final_gene[g_len - 1][2]
        g_name = data.gene_names[i]
        offset = last_offset + last_count
        count = c_idx.shape[0]
        final_gene.append((g_name, offset, count))
    final_exp_np = rfn.unstructured_to_structured(
        np.array(final_exp, dtype=int), np.dtype([('x', np.uint32), ('y', np.uint32), ('count', np.uint16)]))
    genetyp = np.dtype({'names': ['gene', 'offset', 'count'], 'formats': ['S32', np.uint32, np.uint32]})
    final_gene_np = np.array(final_gene, dtype=genetyp)
    h5f = h5py.File(output, "w")
    geneExp = h5f.create_group("geneExp")
    binsz = "bin" + str(data.bin_size)
    bing = geneExp.create_group(binsz)
    geneExp[binsz]["expression"] = final_exp_np  # np.arry([(10,20,2), (20,40,3)], dtype=exptype)
    geneExp[binsz]["gene"] = final_gene_np  # np.arry([("gene1",0,21), ("gene2",21,3)], dtype=genetype)
    if data.attr is not None:
        for key, value in data.attr.items():
            bing["expression"].attrs.create(key, value)
    h5f.attrs.create("version", 2)
    h5f.attrs.create("omics", 'Transcriptomics')
    h5f.close()


write_mid_gef(data, output='./cell_correct_result/cellbin.gef')