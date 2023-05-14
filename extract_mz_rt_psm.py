#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 14:11:52 2022

@author: jwahnzavalet
"""
import pandas as pd
import os
import ntpath
import numpy as np
import re

def convert_to_unmod(match_obj):
    if match_obj.group() is not None:
        return("")
    
def extract_rt_mz (psm_path, top_feature_list=None):
    '''
    Extracts mz and RT for each feature from the psm.tsv file
    ---------------------
    input: path to psm.tsv file, (optional) list of selected features
    returns: matrix of coordinates and list of samples
    '''
    #Read mz from library.tsv file
    coor_col = ['Spectrum', 'Peptide', 'Charge', 'Observed M/Z', 'Retention']
    coor = pd.read_table(psm_path, usecols = coor_col)
    coor['Spectrum'] = [file.split('.')[0] for file in coor['Spectrum'].values]
    coor['Precursor.Id'] = [pep+str(coor.loc[idx, 'Charge']) for idx, pep in enumerate(coor['Peptide']) ]
    coor.drop(columns=['Peptide', 'Charge'], inplace=True)
    coor.set_index('Precursor.Id', inplace=True)
    coor.rename({'Spectrum':'Run', 'Observed M/Z':'mz', 'Retention':'RT'}, axis=1 , inplace =True)
    if top_feature_list != None:
        top_feature_list = [re.sub(r'(\(.*?\))', convert_to_unmod,x) for x in top_feature_list]
        coor = coor.loc[top_feature_list,:]
    
    sample_names = coor['Run'].unique()
    
    return(coor, sample_names)








