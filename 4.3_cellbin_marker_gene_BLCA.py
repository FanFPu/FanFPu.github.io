# -*- coding: utf-8 -*-
'''
Description: 
Author: gongch
Date: 2022-10-25 10:47:44
LastEditTime: 2022-10-25 10:47:45
LastEditors: gongchanghao
E-mail: gongchanghao@genomics.cn
'''
#%%
import os,sys
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from PIL import Image
import random
Image.MAX_IMAGE_PIXELS = None
#%%

def getDefaultColors(n, type = 1):
    if type == 1:
        colors = ["#ff1a1a", "#1B9E77", "#1a1aff", "#ffff1a", "#ff1aff",
                "#ff8d1a", "#7cd5c8", "#c49a3f", "#5d8d9c", "#90353b",
                "#507d41", "#502e71", "#1aff1a", "#c5383c", "#0081d1",
                "#674c2a", "#c8b693", "#aed688", "#f6a97a", "#c6a5cc",
                "#798234", "#6b42c8", "#cf4c8b", "#666666", "#ffd900",
                "#feb308", "#cb7c77", "#68d359", "#6a7dc9", "#c9d73d"]
    elif type == 2:
        if n <= 14:
            colors = ["#437BFE", "#FE6943", "#679966", "#FEC643","#C643FE",
                  "#43D9FE", "#B87A3D", "#43FE69", "#993333", "#7F6699",
                  "#E78AC3", "#333399", "#A6D854", "#E5C494"]
        elif n <= 20:
            colors = ["#87b3d4", "#d5492f", "#6bd155", "#683ec2", "#c9d754",
                  "#d04dc7", "#81d8ae", "#d34a76", "#607d3a", "#6d76cb",
                  "#ce9d3f", "#81357a", "#d3c3a4", "#3c2f5a", "#b96f49",
                  "#4e857e", "#6e282c", "#d293c8", "#393a2a", "#997579"]
        elif n <= 30:
            colors = ["#1aff1a", "#a03259", "#4836be", "#ffff1a", "#ff1aff",
                  "#ff8d1a", "#7cd5c8", "#c49a3f", "#da4f1e", "#90353b",
                  "#507d41", "#502e71", "#1B9E77", "#c5383c", "#0081d1",
                  "#674c2a", "#c8b693", "#aed688", "#f6a97a", "#c6a5cc",
                  "#798234", "#6b42c8", "#cf4c8b", "#666666", "#1a1aff",
                  "#feb308", "#cb7c77", "#68d359", "#8f5b5c", "#c9d73d"]
            #colors = ["#628bac", "#ceda3f", "#7e39c9", "#72d852", "#d849cc",
                  #"#5e8f37", "#5956c8", "#cfa53f", "#392766", "#c7da8b",
                  #"#8d378c", "#68d9a3", "#dd3e34", "#8ed4d5", "#d84787",
                  #"#498770", "#c581d3", "#d27333", "#6680cb", "#83662e",
                  #"#cab7da", "#364627", "#d16263", "#2d384d", "#e0b495",
                  #"#4b272a", "#919071", "#7b3860", "#843028", "#bb7d91"]
        else:
            colors = ["#1aff1a", "#a03259", "#4836be", "#ffff1a", "#ff1aff",
                  "#ff8d1a", "#7cd5c8", "#c49a3f", "#da4f1e", "#90353b",
                  "#507d41", "#502e71", "#1B9E77", "#c5383c", "#0081d1",
                  "#674c2a", "#c8b693", "#aed688", "#f6a97a", "#c6a5cc",
                  "#798234", "#6b42c8", "#cf4c8b", "#666666", "#1a1aff",
                  "#feb308", "#cb7c77", "#68d359", "#8f5b5c", "#c9d73d",
                  "#982f29", "#5ddb53", "#8b35d6", "#a9e047", "#ffd900",
                  "#e0dc33", "#d248d5", "#61a338", "#9765e5", "#69df96",
                  "#7f3095", "#d0d56a", "#371c6b", "#cfa738", "#5066d1",
                  "#e08930", "#83e6d6", "#df4341", "#6a8bd3", "#5d8d9c",
                  "#6ebad4", "#e34c75", "#50975f", "#d548a4", "#badb97",
                  "#b377cf", "#899140", "#564d8b", "#ddb67f", "#292344",
                  "#d0cdb8", "#421b28", "#5eae99", "#ff1a1a", "#406024",
                  "#e598d7", "#343b20", "#bbb5d9", "#975223", "#576e8b",
                  "#d97f5e", "#253e44", "#de959b", "#417265", "#712b5b",
                  "#8c6d30", "#a56c95", "#5f3121", "#8f846e", "#6a7dc9"]
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

def removeBiasGenes(adata):
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
    remove_genes = malat1 | MTgenes | hb_genes | RPgenes | RPgenes2 | CTCgenes | MIRgenes | ACgenes | CTgenes | LINCgenes | ALgenes
    keep = np.invert(remove_genes)
    res = adata[:,keep]
    return res

args = sys.argv
stdata = args[1]
mask = args[2]
outdir = args[3]
resolution = args[4]
os.system("mkdir -p %s"%(outdir))
os.chdir(outdir)
adata = sc.read_h5ad(stdata)

#### Rename cluster
SigGeneral_all_sort = {"Basal_tumor_cell":["KRT5","KRT14","KRT6A"],
                        "muscle_cell":["MYH11","NR2F2","CRYAB","LMOD1","TPPP3"],
                        "lymphoid":["PTPRC","CD247"],
                        "myeloid":["CSF2RA","CSTB","SPP1","APOE","C1QB"],
                        "macrophage":["ADGRE1","CD163", "SLC11A1", "APOC1", "CD86", "CSF1R", "SLCO2B1", "CD68", "F13A1", "CD14", "AIF1", "CD80","FCER1G", "FCGR3A", "TYROBP","LYZ","MS4A7"],
                        "fibroblast":["CALD1","LUM","COL1A1","COL1A2","ACTA2", "SULF1", "CTGF", "TAGLN","TPM1", "GINS1" ,"THY1", "RBP1", "COL1A2", "COL1A1", "C1R", "IGFBP7", "SFRP2", "MGP","C1S", "DCN", "CXCL14", "COL3A1","FAP","COL6A1","PDPN"],
                        "Tcell":["CCL5","CD52","KLRB1","CD2", "CD3D", "CD3E", "CD3G", "CD8A", "CD8B", "GZMK", "CD4","TNFRSF4","IL7R"],
                        "Bcell":["KIF4A","FANCI","CHAF1A", "GTSE1", "ASPM", "SPC25","NCAPG2","POLA2","NCAPD3", "CD19", "CD21", "MS4A1", "CD79A", "CD79B", "BLNK","MZB1"],
                        "DCs":["HLA-DRB1","ITGAX","CD83", "HLA-DMA","HLA-DQB1", "HLA-DPA1","HLA-DPB1","AIF1","LST1","FTL","HLA-C","CD209"],
                        "Epithelial":["HSPA6","S100A2","KRT17","KRT5","SLPI","CXCL17","C1orf56","S100A9","S100A8","IFI44L","EPCAM", "KRT3", "KRT14", "MUC1", "TP63","CDH1","KRT13","KRT17"],
                        "Urothelial_cell":["KRT13","UPK2","UPK1B","UPK3A"],
                        "Epithe_UPK":["GATA3", "PPARG"],
                        "Epithe_CDH12_18":["CDH18","CDH12"],
                        "Epithe_cycling":["TUBA1B","EEF2"],
                        "Epithe_Lum":["HIF1A","CEBPD","BTG1","KRT8","CD9","AQP3","MUC1","KRT18","SLPI","AGR2","LCN2","CLDN4","ANXA1","CD74"],
                        "Epithe_Mam":["PRLR","CLDN4","CSN3","CSN1S1","KRT19","KRT7","KRT8","KRT18"],
                        "myoepithelial":["CD34","CNN1","CNN2","CNN3","IL12B","TP63",'ACTA2','KRT15'],
                        "PVL":["MCAM", "CD146", "ACTA2", "PDGFRB", "LYVE1","RGS5","COL4A1","CALD1","SPARCL1","SPARC","NID1"],
                        "Pericyte":["ACTA2","PDGFRB","MCAM","HIGD1B","ANGPT2","VIM","MFGE8","MYO1B","NOTCH3","COX4I2"],
                        "endothelial":["FLT1","PLVAP","SPARCL1","GNG11","PECAM1", "CD31", "CD34", "HSPG2", "LDB2", "GPR116", "PTPRB", "VWF", "DOCK9", "CDH5", "SELE","VCAM1","ENG"],
                        "endocrine":["CHGB", "CHGA", "TTR", "SCG5", "SLC30A8", "GCG", "CLU","CPE","SCG3","CRYBA2","TM4SF4","SCGN"],
                        "monocyte":["LYZ","MS4A7","CD14"],
                        "tuftCells":["AZGP1", "PLCG2", "HPGDS", "AVIL", "PAEP", "SH2D6", "BMX", "LRMP"],
                        "neutrophils":["A1BG","ALOX5","ASAH1","CD33","CD44","CD63","CTSG","DOCK2","HSPA1B","HSP90AA1"],
                        "mastCells":["RHOH","BTK","FER","GATA2","IL4R","KIT","LCP2","ENPP3","RAC2","LAT2","CD84","LAT","ADGRE2","UNC13D","NDEL1", "CPA3","TPSB2","TPSAB1","MS4A2","SLC18A2","IL1RL1","PTGS1","HPGDS"],
                        "NKcell":["GNLY", "FCER1G", "KLRB1", "KLRC1", "AREG", "XCL1"],
                        "NKT":["NKG7", "GNLY", "GZMA", "GZMB", "FCGR3A", "KLRB1"],
                        "DuctalCell":["CAPS","TPPP3","MIA","RSPH1","PIFO","LCN2","GDF15","AGR3","CETN2","AMBP", "FXYD2"],
                        "plasmocyte":["IGLC2","IGLC3","IGHG3","IGHG1","IGHG4","JCHAIN","IGLC1"],
                        "Cyclegene":["MKI67","CDC20", "CENPF","PTTG1","TOP2A","CCNB1","PCNA"],
                        "astrocyte":["GFAP", "BMPR1B", "CD44", "SLC1A2", "AQP4", "S100B", "GJB7", "ALDH1L1", "ALDOC", "MLC1"],
                        "oligodendrocyte":["MBP", "SOX10", "MOG", "CA2", "CNP", "RTN4", "PLP1", "PLP2", "OPALIN", "OMG","OLIG1", "TNR", "ALCAM", "PLLP"],
                        "MSC":["CD44", "ITGA1","NT5E","THY1"],
                        "mycaf":["RGS5","ACTA2","VIM","CCN2", "COL1A1","COL5A1","COL6A1","TNC","TGFB1","THY1","TAGLN","COL12A1","PDGFRB"],
                        "icaf":["PDGFRA","IL1A","IL1B","IL6","IL11","LIF","CLEC3B","COL14A1","GSN","LY6C1","CXCL12","CXCL14"],
                        "apcaf":["SLPI","SAA3","CD74","H2-Ab1","NKAIN4", "IRF5"],
                        "psc":["DES", "GFAP", "CHRNA1"]}
SigGeneral_select = {}
for i in SigGeneral_all_sort.keys():
    genes=[]
    for j in SigGeneral_all_sort[i]:
        if j in adata.var_names:
            genes.append(j)
    SigGeneral_select[i] = genes

sc.pl.dotplot(adata, SigGeneral_select, groupby=resolution, dendrogram=True, standard_scale='var', color_map="Blues")
plt.savefig("known_marker_genes_dotplot.pdf")
#sc.pl.dotplot(adata, SigGeneral_all_sort, groupby='res.0.8', dendrogram=True, standard_scale='var', color_map="Blues")
#plt.savefig("known_marker_genes_dotplot_res0.8.pdf")
#sc.pl.stacked_violin(adata, SigGeneral_all_sort, groupby='res.0.8', swap_axes=True, dendrogram=True)
#plt.savefig("known_marker_genes_violin.pdf")

adata = removeBiasGenes(adata)
del adata.raw
adata.uns['log1p']["base"] = None
sc.tl.rank_genes_groups(adata, resolution, method='wilcoxon')
plt.rcParams["figure.figsize"] = (7,7)
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
plt.savefig("rank_genes_groups.pdf")

markers = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(3).stack().values.tolist()
markers = list(set(markers))
sc.pl.stacked_violin(adata, markers, groupby=resolution, rotation=90)
plt.savefig("marker_genes_violin.pdf")

degs = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(50)
degs.to_csv("degs.txt",sep = '\t')


sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts'],
             jitter=0.4, multi_panel=True, groupby=resolution)
plt.savefig("violin_genes_count.pdf")

os.system('mkdir -p cluster_split_spatial')
os.chdir("cluster_split_spatial")
flout = open("color.list",'w')
flout.write('1\t#ff0000\n')
flout.write('low_quality\t#ffffff\n')
flout.close()
res = pd.DataFrame(adata.obs, columns = ["x", "y", resolution])
res[resolution]=res[resolution].astype(str)
res.to_csv("bin1clu.txt",sep = '\t',index =False)
clusters = adata.obs[resolution].cat.categories
for i in clusters:
    i = str(i)
    if i == 'low_quality':
        continue
    tmp = res.copy()
    tmp[resolution].loc[tmp[resolution] != i] = 'low_quality'
    tmp[resolution].loc[tmp[resolution] == i] = '1'
    #tmp.loc[tmp["res.0.8"] != i, "res.0.8"] = 'low_quality'
    #tmp.loc[tmp["res.0.8"] == i, "res.0.8"] = '1'
    tmp.to_csv("bin1clu_%s.txt"%(i),sep = '\t',index =False)
    os.system('/jdfsbjcas1/ST_BJ/P21H28400N0232/gongchanghao/bin/cell_bin_plot bin1clu_%s.txt %s color.list cluster_plot_%s.tif'%(i, mask, i))





