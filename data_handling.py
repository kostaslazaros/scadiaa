import os
import requests
import hdf5plugin
from scalex import SCALEX
# from scalex.plot import embedding
# from scalex.metrics import batch_entropy_mixing_score
# from scalex.metrics import silhouette_score
import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
from config import H5AD_PATH, H5AD_CONCATED_PATH, TXT_GZ_PATH, SCREAD_URL, INTEGRATED_PATH
# from matplotlib.pyplot import rc_context


def get_datasets(name: str, prefix: str, postfix: str) -> None:
    if name == '':
        return
    fname = f'{name}{postfix}'
    txt_gz_path = f'./data/scREAD_txt_gz/{fname}'
    if os.path.exists(txt_gz_path):
        print(f'file {fname} already exists in local file system')
        return
    url = f'{SCREAD_URL}{prefix}{fname}'
    response = requests.get(url)
    open(txt_gz_path, 'wb').write(response.content)


def create_h5ad(filename, df_annot) -> None:
    if filename not in df_annot['id'].values:
        return
    h5ad_path = f'{H5AD_PATH}{filename}.h5ad'
    expr_path = f'{TXT_GZ_PATH}{filename}_expr.txt.gz'  
    cell_lbl_path = f'{TXT_GZ_PATH}{filename}_cell_label.txt.gz'
    if os.path.exists(h5ad_path):
        print(f'file {filename} already exists in local file system')
        return
    filename_idx = np.where(df_annot['id'].values == filename)[0][0]
    idv, species, state, region, gender, age, _, _ = df_annot.iloc[filename_idx]
    print(idv, species, state, region, gender, age)
    mtx = sc.read(expr_path, cache=True)
    mtx = mtx.transpose()
    mtx.X = sparse.csr_matrix(mtx.X)
    cell_labels = pd.read_csv(
        cell_lbl_path,
        delimiter='\t',
        names=['cell_name', 'label']
    )
    cell_labels = cell_labels.iloc[1:, :]
    mtx.obs['cell_type'] = list(cell_labels.label.values)
    mtx.obs['species'] = species
    mtx.obs['state'] = state
    mtx.obs['brain_region'] = region
    mtx.obs['gender'] = gender
    mtx.obs['age'] = age
    mtx.obs['data_id'] = filename
    sc.pp.pca(mtx)
    sc.pp.neighbors(mtx)
    sc.tl.umap(mtx)
    writeh5(mtx, h5ad_path)


def merge_data(name_lst: list, merged_h5ad_name):
    anndata_lst = []
    conc_name = f'{H5AD_CONCATED_PATH}{merged_h5ad_name}'
    if os.path.exists(conc_name):
        print(f'File {conc_name} already exists')
        return
    # merged_h5ad_name = f'{list2md5(name_lst)}.h5ad'
    for name in name_lst:
        anndata = sc.read_h5ad(f'{H5AD_PATH}{name}.h5ad')
        anndata_lst.append(anndata)
    merged_anndata = ad.AnnData.concatenate(*anndata_lst, join='inner')
    writeh5(merged_anndata, conc_name)


def writeh5(adata, file_path):
    adata.write_h5ad(
        file_path,
        compression=hdf5plugin.FILTERS['zstd'],
        compression_opts=hdf5plugin.Zstd(clevel=5).filter_options
    )
    

def run_scalex(h5ad_name):
    anndata_corrected = SCALEX(
        f'{H5AD_CONCATED_PATH}{h5ad_name}',
        batch_name='batch',
        min_features=400,
        min_cells=3,
        outdir=f'{INTEGRATED_PATH}{h5ad_name}',
        show=False,
        gpu=7
    )


# def visualize_umap(h5ad_path):
#     anndata = sc.read_h5ad(h5ad_path)
#     with rc_context({'figure.figsize': (6, 6)}):
#         sc.pl.umap(anndata,color=['batch'],legend_fontsize=10)
