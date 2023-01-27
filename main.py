import os
# import hdf5plugin
# import numpy as np
import pandas as pd
# import requests
# import scanpy as sc
import streamlit as st
from st_aggrid import AgGrid, ColumnsAutoSizeMode, GridUpdateMode
from st_aggrid.grid_options_builder import GridOptionsBuilder
# from config import H5AD_PATH, TXT_GZ_PATH, SCREAD_URL
# from scipy import sparse
from data_handling import get_datasets, create_h5ad, merge_data, run_scalex





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
if st.button('Run Integration'):
    st.write(name_lst)
    for elm in name_lst:
        get_datasets(elm, 'expression/', '_expr.txt.gz')
        get_datasets(elm, 'label/', '_cell_label.txt.gz')
        create_h5ad(elm, metadata)
    merged_h5ad_name = merge_data(name_lst)
    st.write('Running SCALEX integration')
    run_scalex(merged_h5ad_name)
    
