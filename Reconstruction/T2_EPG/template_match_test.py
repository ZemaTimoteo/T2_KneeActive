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

from scipy.io import savemat
from tkinter import *
from tkinter import filedialog

sys.path.append('D:/Tiago/Trabalho/2021_2025_PhD/Projects/qMRI_Joint/Code/pythonCode/Reconstruction/T2_EPG')  # Add functions from mytoolbox
from Proc_epg import Proc_epg

def template_match_test(S, X, T2par, B1par):

    # =============================================================================
    # %% --- 1 - normalize the matrices
    # =============================================================================
    X = X / np.linalg.norm( X )

    # X - image, of only n points
    n_points = 3
    aX = X[0:n_points]

    S = abs(S)
    aS = S[0:n_points,:]
    # for i in range(S.shape[1]): # why normalized twice?
    #     S[:,i] = S[:,i] / np.linalg.norm( S[:,i] )
    # =============================================================================
    # %% --- 2 - obtain inner product
    # =============================================================================
    inner_product = X.dot(S)
    ainner_product = aX.dot(aS)

    # =============================================================================
    # %% --- 3 - find the index with the highest inner product for each pixel
    # =============================================================================
    ind_X = np.argmax(abs(inner_product), axis=0)
    aind_X = np.argmax(abs(ainner_product), axis=0)

    # =============================================================================
    # %% --- 4 - Check more common values for dictionary for 1% of max. value
    # =============================================================================
    # # ... 4.1 - Figure ...
    # test_value = 8418
    # origin_X   = T2par[0, ind_X]
    # origin_y   = B1par[0, ind_X]/100
    # test_x     = T2par[0, test_value]
    # test_y     = B1par[0, test_value]/100
    #
    # plt.figure()
    # plt.plot(X)
    # # plt.plot(S[:,ind_X])
    # plt.plot(S[:, test_value])


    # atest_value = 8418
    # aorigin_X   = T2par[0, aind_X]
    # aorigin_y   = B1par[0, aind_X]/100
    # atest_x     = T2par[0, atest_value]
    # atest_y     = B1par[0, atest_value]/100

    # plt.figure()
    # plt.plot(aX)
    # # plt.plot(aS[:,aind_X])
    # plt.plot(aS[:,atest_value])
    # plt.figure()
    # plt.plot(inner_product)


    # # ... 4.1 - Get interval ...
    # aux_interv_Pars = np.where(inner_product > inner_product[ind_X] * .999)
    # interv_Pars     = aux_interv_Pars[0]
    #
    # # ... 4.2 - Get normal parameters ...
    # value   = inner_product[interv_Pars[:]]
    # pars_T2 = T2par[0, interv_Pars[:]]
    # pars_B1 = B1par[0, interv_Pars[:]]/100

    # get t2 values for higher than 0.999 value of inner product.
    # plt.figure()
    # plt.plot(abs(pars_T2))
    # plt.show()

    # aux_x    = T2par[0, ind_X]
    # aux_y    = B1par[0, ind_X]/100
    # aux_value_xy = np.where((pars_T2 == aux_x) & (pars_B1 == aux_y))
    # value_xy     = int(aux_value_xy[0])


    # # ... 4.3 - Figure ...
    # plt.figure()
    # plt.plot(abs(value))
    # labels = pars_T2
    # plt.xticks(labels)
    # plt.plot(value_xy, aux_x, marker="o", markersize=20, markeredgecolor="red",
    #          markerfacecolor="green")
    # plt.show()
    #

    return ind_X