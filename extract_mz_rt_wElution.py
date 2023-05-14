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

def extract_rt_mz (diann_path, library_path, top_feature_list):
    '''
    Extracts the m/z and RT elution profile for each feature
    '''
    #Read RT from diann-output file
    RT_col = ['RT.Start', 'RT.Stop', 'Precursor.Id', 'Run']
    RT = pd.read_table(diann_path, usecols = RT_col)
    RT.set_index('Precursor.Id', inplace=True)
    RT = RT.loc[top_feature_list,:]
    
    #Read mz from library.tsv file
    mz_col = ['PrecursorMz', 'ModifiedPeptideSequence', 'PrecursorCharge']
    mz_raw = pd.read_table(library_path, usecols = mz_col)
    mz = pd.DataFrame({'Precursor.Id' : mz_raw['ModifiedPeptideSequence'] + mz_raw['PrecursorCharge'].astype('string'),
                      'PrecursorMz' : mz_raw['PrecursorMz']})
    mz = mz.drop_duplicates(keep='first')
    mz.set_index('Precursor.Id', inplace=True)
    mz = mz.loc[top_feature_list,:]
    
    sample_names = RT['Run'].unique()
    
    uni_pep = pd.DataFrame({'idx':mz.index, 
                            'mz' : mz['PrecursorMz']})
    uni_pep.set_index('idx', inplace=True)
    
    dic = {}
    
    #Extract rt/mz of each peptide for each run/sample
    for sample in sample_names:
        sample_df = RT.iloc[np.where(RT['Run']==sample)][['RT.Start', 'RT.Stop']]
        
        RT_sample = []
        for pep in uni_pep.index:
            if pep in sample_df.index:
                RT_sample.append(sample_df.loc[pep,['RT.Start', 'RT.Stop']].tolist())
            
            else:
                RT_sample.append([np.nan,np.nan])
        
        dic['RT_'+sample] = RT_sample
        print(dic)
    
    RT_samples_df = pd.DataFrame(dic)
    RT_samples_df.index = uni_pep.index.values
    uni_pep= pd.concat([uni_pep, RT_samples_df], axis=1)
    
    return(uni_pep)
