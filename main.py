import os
import hdf5plugin
import numpy as np
import pandas as pd
import requests
import scanpy as sc
import streamlit as st
from st_aggrid import AgGrid, ColumnsAutoSizeMode, GridUpdateMode
from st_aggrid.grid_options_builder import GridOptionsBuilder
from config import H5AD_PATH, TXT_GZ_PATH, SCREAD_URL
from scipy import sparse

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

    # Start from here
    mtx = sc.read(expr_path, cache=True)
    mtx = mtx.transpose()
    mtx.X = sparse.csr_matrix(mtx.X)

    cell_labels = pd.read_csv(
        cell_lbl_path,
        delimiter='\t',
        names=['cell_name', 'label']
    )
    cell_labels = cell_labels.iloc[1:, :]

    mtx.obs['cell type'] = list(cell_labels.label.values)
    mtx.obs['species'] = species
    mtx.obs['state'] = state
    mtx.obs['brain region'] = region
    mtx.obs['gender'] = gender
    mtx.obs['age'] = age

    mtx.write_h5ad(
        h5ad_path,
        compression=hdf5plugin.FILTERS['zstd'],
        compression_opts=hdf5plugin.Zstd(clevel=5).filter_options
    )


metadata = pd.read_csv('scREAD_meta.csv')


gd = GridOptionsBuilder.from_dataframe(metadata)
gd.configure_selection(selection_mode='multiple', use_checkbox=True)
gridoptions = gd.build()

grid_table = AgGrid(
    metadata,
    height=350,
    gridOptions=gridoptions,
    update_mode=GridUpdateMode.SELECTION_CHANGED,
    columns_auto_size_mode=ColumnsAutoSizeMode.FIT_CONTENTS
)

# st.write('## Selected')
selected_rows = grid_table["selected_rows"]
# st.dataframe(selected_rows)
if st.button('Run Integration'):
    st.write([row['id'] for row in selected_rows])
    for elm in [row['id'] for row in selected_rows]:
        get_datasets(elm, 'expression/', '_expr.txt.gz')
        get_datasets(elm, 'label/', '_cell_label.txt.gz')
        create_h5ad(elm, metadata)
