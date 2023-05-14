#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 22 12:49:52 2022

@author: Johanna Wahn
"""

import os
import os.path
import pandas as pd
import numpy as np
import pycombat
from sklearn.pipeline import Pipeline
from feature_preprocessing import MedianNormalizeIntensityMatrix, FilterPeptideMissingness, FillNanWithMin
from DIANN_to_TRIC import FormatDIANN, TricToPepDf


class DIANN_to_expr_matrix:
    """Processes Fragpipe DIANN output matrix into a peptide expression matrix
        
    Returns:
        -Normalized, batch corrected peptide expression matrix
        
    """
    def __init__(self):
        self = self
    
    def fit(self, X, y=None):
        return self

    def transform(self, data_path, meta_path, save_file_name=None, batch_var='location', batch_corr='True'):
        X = self.load_data_as_TRIC(data_path)
        meta = self.load_meta(meta_path)
        X = self.merge_to_meta(X, meta)
        X = self.get_expr_matrix_normalized(X)
        if batch_corr:
        	X = self.batch_correction(X, batch_var, meta)
        if save_file_name!= None: 
            X.to_csv(save_file_name)    
        else: 
              return X 
    def load_data_as_TRIC(self, data_path):
        """Loads and reformats Fragpipe Matrix to TRIC-like one
        
        This function is modeled as a scikit transformer
        
        Args:
            - Filepath to the fragpipe output file to be processed
            
        Returns: TRIC formatted Fragpipe matrix
            
        """
        load_DIANN = FormatDIANN(intensity_column="Precursor.Quantity")
        qp_long =load_DIANN.transform(filename=data_path)
        
        sample_names = pd.DataFrame(qp_long['filename'].unique())
    
        print('There are {} samples in the quantitative matrix'.format(sample_names.shape[0]))
        print('Samples appearing more than once:')
        print(sample_names[sample_names.duplicated()])
        return(qp_long)
    
    def load_meta(self,meta_path):
        """Loads and reformats metafile 
        
        Args:
            - Filepath to the metafile
            
        Returns: Reformatted metafile
            
        """
        meta = pd.read_csv(meta_path)
        
        def extract_run_order(runname):
            return (runname[runname.find('_')+2:])
        
        def add_locations(label):
            if "Basel" in label:
                return ("Basel")
            return ("Zurich")
        
        meta['run_order'] = meta['File.Name'].apply(lambda x: extract_run_order(x))
        meta['location'] = meta['Run.Label'].apply(lambda x: add_locations(x))
        meta = meta.rename(columns={'Condition':'class',
                                    'Run.Label': 'runlabel',
                                    'File.Name': 'filename'})
        meta.sort_values('run_order', inplace = True)
        
        meta.index = meta['filename'].values
        meta.sort_values('run_order', inplace = True)
        meta.reset_index(drop=True, inplace=True)
        meta['rank'] = meta.index
        return(meta)
        
    
    def merge_to_meta(self, X, meta):
        """Merges metafile data to main matrix 
        
        Args:
            - Peptide matrix
            - Metafile
            
        Returns: Reformatted peptide matrix
            
        """
        qp_l = X.merge(meta)
        return(qp_l)
    
    def get_expr_matrix_normalized(self, X):
        """Transorms Fragpipe matrix to normalized peptide expression matrix 
        
        Args:
            - Peptide matrix
        
        Returns: Normalized peptide expression matrix 
            
        """
        pep_mtx_t = TricToPepDf().fit_transform(X)
        
        pep_mtx_t.replace(-np.inf, 0, inplace=True)
        preprocessing_pipe = Pipeline([('filter_peptide_missigness', FilterPeptideMissingness(0.5)),
                                       ('filter_sample_missigness', FilterPeptideMissingness(0.5)),                               
                                       ('median_normalize', MedianNormalizeIntensityMatrix()),
                                       ('fill_nan_with_min', FillNanWithMin(randomize=True))
                                      ])
    
        pep_mtx_t_norm = preprocessing_pipe.fit_transform(pep_mtx_t)
        pep_mtx_norm = pep_mtx_t_norm.transpose()
        return(pep_mtx_norm)
    
    def batch_correction (self, data, batch, meta):
        """Batch correction for a specified variable (set by default to location)
        
        Args:
            - Peptide expression matrix
            - Batch variable to be corrected 
            - Metafile
        
        Returns: Batch corrected peptide expression matrix 
            
        """
        X = data.copy()
        factors = meta.set_index('filename')
        X = X.join(factors)
        y = pd.DataFrame(X[batch])
        X = X.drop(factors.columns, axis=1)
        
        combat = pycombat.Combat()
        
        X_cbt = pd.DataFrame(combat.fit_transform(X.values, y))
        X_cbt.columns = X.columns
        X_cbt.columns.name = 'FullPeptideName'
        X_cbt.index = X.index
        return(X_cbt)

    


