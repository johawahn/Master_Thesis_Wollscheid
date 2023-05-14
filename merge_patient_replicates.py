#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 15:47:37 2022

@author: jwahnzavalet
"""
import pandas as pd
import numpy as np


def reduce_expr_matrix(table, patients):
    '''
    Function that takes the peptide/protein expression matrix and averages the measurements of 
    files corresponding to replicates of the same patient
    
    Parameters
    ----------
    table : table of peptide/protein expression with file names as rows
    patients : patient annotation file containing patient IDs and their respective files

    Returns
    -------
    Expression table with merged patient replicates.

    '''
    def filter_patients(used_files, meta_files):
        '''
        Function that reduces the meta file to only the files present in the table of measurements
        
        Parameters
        ----------
        used_files : list of filenames of interest
        
        meta_files : meta file of patient annotations

        Returns
        -------
        Meta files of only the files of interest

        '''
        f_meta_files = meta_files.set_index("File")
        f_meta_files = f_meta_files.loc[used_files]
        
        return(f_meta_files.reset_index())
    
    patients = filter_patients(table["filename"], patients)
    
    patients = patients.set_index("PatientID")
    table = table.set_index('filename')
    
    data_to_append = []
    id_to_append =[]
    for id in np.unique(patients.index.values):
        files = patients.loc[id].values
        id_to_append.append('Patient'+str(id))
        if len(files) > 1:
            data_to_append.append(table.loc[patients.loc[id,'File']].mean(axis=0).values)
    
        else:
            data_to_append.append(table.loc[patients.loc[id,'File']].values)
    
    reduced_t = pd.DataFrame(data_to_append)
    reduced_t.columns = table.columns
    reduced_t.index = id_to_append
    
    return(reduced_t)


