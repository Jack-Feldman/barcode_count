#!/usr/bin/env python3
"""
Functions used in barcode_extract_widget.py script

"""

import pandas as pd
import numpy as np

def get_cell_replicates(data: pd.DataFrame) -> dict:
    """
    Returns dictionary
    {cell type : [list of replicate names of that cell type]}
    """
    types = {}
    for name in data.columns:
        cell_type = name.split(".")[0]
        if cell_type not in types:
            types[cell_type] = [name]
        else:
            types[cell_type].append(name)

    return types

def calc_metrics(data: str):
    """
    Input
        File path to normcounts.csv file
    Output
        New df with calculated std and average metrics
    """
    data = pd.read_csv(data, index_col=[0])
    # Create new df
    cell_data = pd.DataFrame(index = data.index)
    compiled_data = pd.DataFrame(index = data.index)
    cell_types = get_cell_replicates(data)

    for cell_type in cell_types:

        cell_mean = data[cell_types[cell_type]].mean(axis=1)
        cell_sd = data[cell_types[cell_type]].std(axis=1)

        for output_file in [cell_data, compiled_data]:
            output_file[cell_type + '_avg'] = cell_mean
            output_file[cell_type + '_std'] = cell_sd

        #Add norm counts into compiled data file
        compiled_data[cell_types[cell_type]] = data[cell_types[cell_type]]

    return cell_data, compiled_data

def return_pcr_bias_df() -> pd.DataFrame:
    """
    Returns multi-indexed df for counting pcr bias like below:

        		A	C	G	T
region	position
region_1	1	2.0	1.0	1.0	1.0
            2	1.0	2.0	1.0	1.0
            3	2.0	1.0	1.0	1.0
            4	1.0	2.0	1.0	1.0
region_2	1	1.0	1.0	2.0	1.0
            2	1.0	1.0	1.0	2.0
            3	1.0	2.0	1.0	1.0
            4	2.0	1.0	1.0	1.0
region_3	1	1.0	1.0	1.0	2.0
            2	1.0	1.0	2.0	1.0
            3	1.0	2.0	1.0	1.0
    """

    arrays = [['region_1']*4+['region_2']*4+['region_3']*3,
              [1,2,3,4]*2+[1,2,3]]
    index = list(zip(*arrays))
    index = pd.MultiIndex.from_tuples(index, names = ['region', 'position'])

    return pd.DataFrame(data = np.ones((11,4)), index=index, columns=['A', 'C', 'G', 'T'])
