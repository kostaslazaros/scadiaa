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
st.set_option('deprecation.showPyplotGlobalUse', False)


def visualize_umap(h5ad_path, color_key):
    anndata = sc.read_h5ad(h5ad_path)
    # umap1 = anndata.obsm['X_umap'][:,0]
    # umap2 = anndata.obsm['X_umap'][:,1]
    fig, ax = plt.subplots()
    
    # ax.scatter(umap1, umap2)
    # st.pyplot(fig)
    st.pyplot(sc.pl.umap(anndata,color=[color_key],legend_fontsize=10))


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

# st.radio(
#     "Set selectbox label visibility 👉",
#     key="umap_key",
#     options=["data_id", "state", "cell_type", "gender"],
#     )
co1, co2, co3, co4 = st.columns(4)
with co1:
    data_id = st.checkbox('data_id')
with co2:
    state = st.checkbox('state')
with co3:
    cell_type = st.checkbox('cell_type')
with co4:
    gender = st.checkbox('gender')
# st.write(st.session_state.umap_key)
if st.button('Run Integration'):

    st.write(name_lst)
    for elm in name_lst:
        get_datasets(elm, 'expression/', '_expr.txt.gz')
        get_datasets(elm, 'label/', '_cell_label.txt.gz')
        create_h5ad(elm, metadata)
    merged_h5ad_name = merge_data(name_lst)
    st.write('Running SCALEX integration')
    run_scalex(merged_h5ad_name)
    if data_id:
        visualize_umap('./A58d3c8988e89f8c/adata.h5ad', "data_id")
    if state:
        visualize_umap('./A58d3c8988e89f8c/adata.h5ad', "state")
    if cell_type:
        visualize_umap('./A58d3c8988e89f8c/adata.h5ad', "cell_type")
    if gender:
        visualize_umap('./A58d3c8988e89f8c/adata.h5ad', "gender")
    
  