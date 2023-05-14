#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 16:15:15 2023

@author: jwahnzavalet
"""

import os 
import argparse
from fasta2csv import prosit_fasta2csv

def csv_to_prosit (server_nbr):
    os.system('env ALL_PROXY=socks5h://localhost:9999 curl -F "peptides=@peptides.csv" ' 
                           + server_nbr +'predict/generic > prosit_lib.csv')
    return()

def prosit_to_mzML (synth_config):
    os.system('synthedia --config ' + synth_config)
    return()

def parse_args():
    """

    Returns
    -------
    args : TYPE
        Default argument parser for console execution.

    """
    parser = argparse.ArgumentParser(description='Convert a fasta file to mzML synthetic data')
    parser.add_argument('--ff', '--fasta_file', type=str, default = os.path.join(os.getcwd(),'test_file.fasta'), help='path to the faster file: "path/to/file', required=True)
    parser.add_argument('--ce', '--collision_energy', type=int, default = 28, help='integer between 10 to 50 identifying keV used for fragmentation')
    parser.add_argument('--cs', '--charge_states', type=list, default = [2], help='list of charge states from which randomly will be selected')
    parser.add_argument('--ps', '--prosit_server', type=str, help='string of ruler server id where prosit server is running', required=True)
    parser.add_argument('--sc', '--synthedia_configurations', type=str, help='path to synthedia configuration file', required=True)
    parser.add_argument('--sim_all', '--simulate_all_charges', type=bool, default='False', help='True or False to simulate all charges for each peptide')
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parse_args()
    prosit_fasta2csv(args.ff, args.ce, args.cs, args.sim_all) 
    print('================================================================================\nFASTA to csv conversion done\n================================================================================')
    csv_to_prosit(args.ps)
    print('===================================================================================\nProsit spectral library done\n================================================================================')
    prosit_to_mzML(args.sc)
    print("================================================================================\nmzML Synthedia file done\n================================================================================")
