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
    
        
    #Read mz from library.tsv file
    coor_col = ['Spectrum', 'Peptide', 'Charge', 'Observed M/Z', 'Retention', 'Entry Name']
    coor = pd.read_table(psm_path, usecols = coor_col)
    coor['Spectrum'] = [file.split('.')[0] for file in coor['Spectrum'].values]
    coor['Precursor.Id'] = [pep+str(coor.loc[idx, 'Charge']) for idx, pep in enumerate(coor['Peptide']) ]
    coor.drop(columns=['Peptide', 'Charge'], inplace=True)
    coor.set_index('Precursor.Id', inplace=True)
    coor.rename({'Spectrum':'Run', 'Observed M/Z':'mz', 'Retention':'RT', 'Entry Name':'ProteinGroup'}, axis=1 , inplace =True)
    if top_feature_list != None:
        top_feature_list = [re.sub(r'(\(.*?\))', convert_to_unmod,x) for x in top_feature_list]
        coor = coor.loc[top_feature_list,:]
    
    sample_names = coor['Run'].unique()
    
    return(coor, sample_names)

#feature_list = pd.read_csv('/Volumes/Untitled/230112_RFE_REDO_221219_Fragpipe/230112_peptide_expr_matrix_w_classes_most_contr_feat.csv' , header=None, usecols=[0])[0].to_list()

#path = '/Users/johannawahn/Desktop/Wollscheid_Lab/psm.tsv'
#coordinates, files = extract_rt_mz(path)

#coordinates.to_csv('/Users/johannawahn/Desktop/230215_interactive_plot/coordinates.csv')







