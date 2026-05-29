# -*- coding: utf-8 -*-
"""
Created on Wed Ago 09 13:37:25 2021
Definition:
Get the template match between the dictionary and the data obtain from the scanner

Functions used:
-

Inputs:
- S: Dict_knee_shortTE_norm - Dictionary normalized for MRI scanner
- X: Image data from MRI Scanner

Outputs:
- ind_X: return the index's of T2 that optimize the data with the dictionary

@author: tfernandes

Needs:
 -
"""

# =============================================================================
#%% --- 0 - Import functions
# =============================================================================

import os
import scipy.io
import matplotlib
import tkinter
import numpy as np
matplotlib.use('TKAgg')

import matplotlib.pyplot as plt
import sys
import math
import scipy.linalg as sc

from scipy.io import savemat
from tkinter import *
from tkinter import filedialog
from matplotlib.colors import ListedColormap

sys.path.append('D:/Tiago/Trabalho/2021_2025_PhD/Projects/qMRI_Joint/Code/pythonCode/Reconstruction/T2_EPG')  # Add functions from mytoolbox
from Proc_epg import Proc_epg
from dict_pars_generator_seq import dict_pars_generator_seq

plt.switch_backend('agg')

def template_match_vComplex(S, X, T2par, B1par):

    # =============================================================================
    # %% --- 1 - normalize the matrices
    # =============================================================================
    X = X / np.linalg.norm( X )

    # X - image, of only n points
    n_points = 3
    aX = abs(X[0:n_points])

    aS = S[0:n_points,:]
    auxS = np.zeros((S.shape), dtype=np.complex64 )
    for i in range(S.shape[1]): # why normalized twice?
        auxS[:,i] = S[:,i] / np.linalg.norm( S[:,i] )
    # =============================================================================
    # %% --- 2 - obtain inner product
    # =============================================================================
    inner_product = abs(X).dot(abs(S))

    # aux_inner_product = X.dot(auxS)
    # ainner_product = aX.dot(aS)

    # =============================================================================
    # %% --- 3 - find the index with the highest inner product for each pixel
    # =============================================================================

    if len(inner_product) == 0:
        a=2
    ind_X     = np.argmax((inner_product), axis=0)

    # =============================================================================
    # %% --- 4 - Plot Tests
    # =============================================================================
    #
    # plt.figure()
    # plt.title('Dictionary')
    #plt.ylabel('a.u.')
    #plt.xlabel('#Echos')
    #
    # plt.figure()
    # Y = np.linspace(0,X.shape[0]-1,X.shape[0])
    # plt.plot(np.squeeze(abs(S[:,ind_X])))
    # plt.scatter(Y,abs(X))
    # plt.title('Signal with Dictionary Best Fit | T2='+str(abs(T2par[:,ind_X][0])) + 'ms & B1=' +str(abs(B1par[:,ind_X][0]))+'%')
    # plt.ylabel('a.u.')
    # plt.xlabel('#Echos')
    # #
    # np.argwhere((T2par[:, :][0] == 30) & (B1par[:, :][0] == 114))
    # #
    # plt.figure()
    # Y = np.linspace(0,X.shape[0]-1,X.shape[0])
    # plt.scatter(Y,X)
    # plt.plot(S[:,ind_X],label='Best Fit | B1='+str(B1par[:,ind_X][0]) + '%')
    # new_idx = np.argwhere((T2par[:, :][0] == 30) & (B1par[:, :][0] == 94))
    # plt.plot(S[:,new_idx[0][0]],label='Fit for | B1='+str(B1par[:,new_idx[0][0]][0]) + '%')
    # plt.title('Signal with Dictionary Fit | T2='+str(T2par[:,ind_X][0]) + 'ms')
    # plt.ylabel('a.u.')
    # plt.xlabel('#Echos')
    # plt.legend(loc="upper right")

    return ind_X