import os
# import hdf5plugin
# import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# import requests
# import scanpy as sc
import streamlit as st
import scanpy as sc
from st_aggrid import AgGrid, ColumnsAutoSizeMode, GridUpdateMode
from st_aggrid.grid_options_builder import GridOptionsBuilder
# from config import H5AD_PATH, TXT_GZ_PATH, SCREAD_URL
# from scipy import sparse
from data_handling import get_datasets, create_h5ad, merge_data, run_scalex
from hash_encode import list2md5
from config import INTEGRATED_PATH, H5AD_PATH
st.set_page_config(layout='wide')
st.set_option('deprecation.showPyplotGlobalUse', False)


def visualize_umap(h5ad_path, color_key):
    anndata = sc.read_h5ad(h5ad_path)
    # umap1 = anndata.obsm['X_umap'][:,0]
    # umap2 = anndata.obsm['X_umap'][:,1]
    fig, ax = plt.subplots()
    fig = sc.pl.umap(anndata, color=[color_key], legend_fontsize=10)
    
    # ax.scatter(umap1, umap2)
    st.pyplot(fig)
    # st.pyplot(sc.pl.umap(anndata,color=[color_key],legend_fontsize=10))


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
name_lst = [row['id'] for row in selected_rows]

co1, co2, co3, co4, co5 = st.columns(5)
with co1:
    data_id = st.checkbox('data_id')
with co2:
    state = st.checkbox('state')
with co3:
    cell_type = st.checkbox('cell_type')
with co4:
    gender = st.checkbox('gender')
with co5:
    species = st.checkbox('species')

if st.button('Run Integration') and len(name_lst) > 0:

    st.write(name_lst)
    for elm in name_lst:
        get_datasets(elm, 'expression/', '_expr.txt.gz')
        get_datasets(elm, 'label/', '_cell_label.txt.gz')
        create_h5ad(elm, metadata)
    if len(name_lst) > 1:
        merged_h5ad_name = f'{list2md5(name_lst)}.h5ad'
        merge_data(name_lst, merged_h5ad_name)
        integrated_path = f'{INTEGRATED_PATH}{merged_h5ad_name}/adata.h5ad'
        if not os.path.exists(integrated_path):
            st.write('Running SCALEX integration')
            run_scalex(merged_h5ad_name)
    else:
        integrated_path = f'{H5AD_PATH}{name_lst[0]}.h5ad'

    if data_id:
        visualize_umap(integrated_path, "data_id")
    if state:
        visualize_umap(integrated_path, "state")
    if cell_type:
        visualize_umap(integrated_path, "cell_type")
    if gender:
        visualize_umap(integrated_path, "gender")
    if species:
        visualize_umap(integrated_path, "species")
