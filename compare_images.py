#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 17 10:15:10 2023

@author: jwahnzavalet
"""

from PIL import Image, ImageOps
from collections import Counter
import numpy as np
import os 
import multiprocessing as mp
import time
from numpy import savetxt

def L2Norm(H1,H2):
    '''
    Calculates the euclidean distance between the images.

    Parameters
    ----------
    H1 : Flat vector of image 1
    H2 : Flat vector of image 2

    Returns
    -------
    Euclidean distance

    '''
    distance = sum([np.square(x1-x2) for x1,x2 in zip(H1, H2)])
    return np.sqrt(distance)

def count_hist_vect(flat_array):
    '''
    Generates the Count-Histogram-Vector

    Parameters
    ----------
    flat_array : 1D array of image

    Returns
    -------
    Counter of values 

    '''
    RH = Counter(flat_array)
    H = [RH[i] if i in RH.keys() else 0 for i in range(256)]

    return(H)

def compare_2Images(img1, img2):
    '''
    

    Parameters
    ----------
    img1 : First image to be compared
    img2 : Second image to be compared

    Returns
    -------
    Distance value between the two images

    '''
    shape = (1193,3000)
    
    img1 = np.resize(np.asarray(ImageOps.grayscale(Image.open(img1))),shape)
    img2 = np.resize(np.asarray(ImageOps.grayscale(Image.open(img2))),shape)
    
    flat_array_1 = img1.flatten()
    flat_array_2 = img2.flatten()

    H1 = count_hist_vect(flat_array_1)
    H2 = count_hist_vect(flat_array_2)
    
    return(L2Norm(H1,H2))

def create_dist_matrix_parallel(path, chunks_size=5):
    '''
    Function to test similarity between all images fo a given folder. 
    This function parallelizes the calculations and takes as many cores that there
    are available.

    Parameters
    ----------
    path : path to folder with all images

    Returns
    -------
    dist_matrix : euclidean distance matrix of all images

    '''
    os.chdir(path)
    files_list = os.listdir(path)
    
    n = len(files_list)
    dist_matrix = np.zeros((n,n))

    # define the number of processes to use
    num_processes = mp.cpu_count()

    # create a pool of processes
    pool = mp.Pool(num_processes)

    # split the range into chunks
    chunks = [(i, min(i+chunks_size, n)) for i in range(0, n, chunks_size)]

    # process each chunk in parallel
    for i, chunk in enumerate(chunks):
        results = pool.starmap(compare_2Images, [(files_list[j], files_list[k]) for j in range(chunk[0], chunk[1]) for k in range(n)])
        results_as_matrix = np.array(results).reshape(chunk[1]-chunk[0],n)
        dist_matrix[chunk[0]:chunk[1],:] = results_as_matrix
    

    pool.close()
    pool.join()
    
    return dist_matrix


def create_dist_matrix(file_list):
    '''
    Function to test similarity between all images fo a given folder. 
    This function does not parallelize, so its adviced to use the function create_dist_matrix_parallel
    if the job contains a large amount of images.
    
    Parameters
    ----------
    path : path to folder with all images

    Returns
    -------
    dist_matrix : euclidean distance matrix of all images

    '''
    n = len(file_list)
    dist_matrix = np.zeros((n,n))    # initialize distance matrix to a square of zeros
    for i in range(n):
        for j in range(i, n):
            dist_matrix[i,j] = compare_2Images(file_list[i], file_list[j])
            dist_matrix[j,i] = dist_matrix[i,j] 
    
    return(dist_matrix)

def replicates_test(dist_matrix):
    '''
    Indicates if there are images that are identical

    Parameters
    ----------
    dist_matrix : Euclidean distance matrix

    Returns
    -------
    None.

    '''
    zeros = np.where(dist_matrix==0)
    
    for idx in range(len(zeros)):
        if zeros[0][idx]!=zeros[1][idx]:
            print(zeros[0][idx],zeros[1][idx])
            
            

