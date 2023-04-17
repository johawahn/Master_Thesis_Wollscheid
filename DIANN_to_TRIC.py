#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 12:36:02 2022

@author: Johanna Wahn and Patrick Pedrioli 
"""
"""Utility functions for pre-processing DIA-NN quantitative proteomic matrices
"""

import pandas as pd
import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin


def quantity_normalised_selection():
    """Selects the Intensity normalised of not depending on user input
    

    Return:
        Name of the DIANN column with normalised or not normalised intensity
    """
    
    normalised = input("Do you want the intensity normalised? [y/n]")
    global norm_int
    is_input = True
    try:
        normalised in ['y', 'n']
    except ValueError:
        is_input = False
        
    if (is_input & (normalised == 'y')):
        norm_int = 'Precursor.Normalised'
        return(norm_int)
    
    elif (is_input & (normalised == 'n')):
        norm_int = 'Precursor.Quantity'
        return(norm_int)
    
    else:
        print("Invalid input!")       
    



class FormatDIANN:
    """Reformat Fragpipe Matrix to TRIC-like one
    
    This function is modeled as a scikit transformer
    
    Args:
        - X: A Fragpipe quantitative peptide matrix.
        Requires;
        - Precursor.Id
        - Protein.Ids
        - Intensity Columns (Precursor Normalised)
        
    """
    def __init__(self, intensity_column):
        self.intensity_column = intensity_column
    
    
    def fit(self, X, y=None):
        return self
    
    def transform(self, filename):
      X = self.load_DIANN_matrix(filename)
      X = self.rename_columns(X)
      return X
    
    def load_DIANN_matrix(self, filename):
        """Transforms a DIANN file into a pandas DataFrame
        
        Args:
            Filename (string): Full path to the DIA-NN matrix
            
        Return: 
            A quantitative feature matrix loaded from DIA-NN output file
        """
        
        feature_matrix = pd.read_csv(filename,
                                     sep="\t",
                                     usecols=[
                                      "Run",
                                      "Protein.Ids",
                                      "Precursor.Id",
                                      "Stripped.Sequence",
                                      "Precursor.Charge",
                                      "Q.Value",
                                      self.intensity_column,
                                      "RT",
                                      ])
        return feature_matrix
    
    def rename_columns(self, X):
        X = X.rename(columns={"Run": "filename", 
                              "Protein.Ids": "ProteinName",
                              "Precursor.Id": "FullPeptideName",
                              "Stripped.Sequence": "Sequence",
                              "Precursor.Charge": "Charge",
                              "Q.Value": "Qval",
                              self.intensity_column: "Intensity"})
        return(X)


class TricToPepDf(BaseEstimator, TransformerMixin):
    """Converts TRIC feature alignment matrix to a peptide level matrix.
    Integrates multiple charge states for a given peptide to generate a
    peptide x sample matrix of intensities.
    
    This function is modeled as a scikit transformer
    
    Args:
      log_data: wether intensities should be logged during conversion
      X: The TRIC feature alignment matrix.
    
    Return:
      A peptide level intensity matrix (peptide x sample)
    """
    def __init__(self, log_data=True):
      self.log_data = log_data
    
    def fit(self, X, y=None):
      return self
    
    def transform(self, X):
      pep_df = (
        X.groupby(["filename", "ProteinName", "FullPeptideName", "Sequence"])
        .Intensity.sum()
        .reset_index()
      )
      if self.log_data:
        pep_df['Intensity'] = np.log10(pep_df['Intensity'])
      X = pep_df.pivot(index='FullPeptideName', columns='filename', values='Intensity')
      return X 


    

    
    
        
