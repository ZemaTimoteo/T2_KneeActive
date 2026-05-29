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


def template_match(S,X):

    # =============================================================================
    # %% --- 1 - normalize the matrices
    # =============================================================================
    X = X / np.linalg.norm( X )
    S = abs(S)

    # for i in range(S.shape[1]): # why normalized twice?
    #     S[:,i] = S[:,i] / np.linalg.norm( S[:,i] )
    # =============================================================================
    # %% --- 2 - obtain inner product
    # =============================================================================
    inner_product = X.dot(S)
    # 1/ std (vial para cada echo) -> normalizar pelo sumatório dos pesos que tem de ser igual a 1.

    # =============================================================================
    # %% --- 3 - find the index with the highest inner product for each pixel
    # =============================================================================
    ind_X = np.argmax(abs(inner_product), axis=0)


    # # ... 4.3 - Figure ...
    # plt.figure()
    # plt.plot(X)
    # plt.plot(S[:,ind_X])
    #
    # plt.show()
    #
    #
    # plt.figure()
    # plt.plot(inner_product)

    # return ind_X, aind_X
    return ind_X