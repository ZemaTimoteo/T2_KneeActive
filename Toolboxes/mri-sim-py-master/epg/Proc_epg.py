#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Ago 2 14:10:22 2021
Procedure for EPG simulations
@author: tfernandes

Needs:
 - Add toolbox from https://github.com/ut - mri-sim-py in a toolbox folder
 -
"""

# =============================================================================
#%% --- 0 - Import functions
# =============================================================================

PC = 0 # PC = 1 - Seia OR PC = 0 - My PC

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
from roipoly import RoiPoly


if PC == 0:
    sys.path.append('D:/Tiago/Trabalho/2021_2025_PhD/Projects/qMRI_Joint/Code/pythonCode/Toolboxes/mri-sim-py-master/epg') # Add toolboxes from https://github.com/utcsilab
    from epgcpmg import FSE_TE
    sys.path.append('D:/Tiago/Trabalho/2021_2025_PhD/Projects/qMRI_Joint/Code/pythonCode/Reconstruction/T2_EPG')  # Add functions from mytoolbox
    from epg_cpmg_LsB import epg_cpmg_LsB

def Proc_epg(exc_pulse, refoc_pulse, exc_phase, refoc_phase, T1, T2, dTE, ETL):
    """
    Input:
        :param exc_pulse:       Magnitude / Profile of excitatory pulse
        :param refoc_pulse:     Magnitude / Profile of Refocusing pulse
        :param exc_phase:       Phase of Excitatory Pulse
        :param refoc_phase:     Phase of Refocusing Pulse
        :param T1:              T1 value in s
        :param T2:              T2 value in s
        :param dTE:             Echo Spacing in s
        :param ETL:             Echo Train Lenght

    Output:
        :return: Dictionary simulation Parameters
    """

    # =============================================================================
    #%% --- 1 - Set # parameters
    # =============================================================================

    echoes = np.zeros((ETL,exc_pulse.shape[0]), dtype=float) # % zero matrix ETL * Length exc_pulse
    soma   = np.zeros( ( ETL, int((T2.shape)[0]),  exc_phase.shape[0] ), dtype=float)

    # =============================================================================
    # %% --- 2 - Generate epg dictionary using epg_cpmg.m function
    # =============================================================================
    # % -- Get value for dic. from CPMG condition
    for j in range( exc_phase.shape[0] ):
        for i in range( int((T2.shape)[0]) ):
            for z in range(  exc_pulse.shape[0] ): # 1 -->numero de TR(128)
                print('   ->  T2_i:', i+1, '/', int((T2.shape)[0]), ' tests')
                # s = epg_cpmg_LsB(exc_pulse[z], exc_phase[j], refoc_pulse[z,:], ETL, T1, T2[i], dTE, 0, refoc_phase)
                s = epg_cpmg_LsB(exc_pulse[z], exc_phase[j], refoc_pulse[z], ETL, T1, T2[i], dTE, 0, refoc_phase)
                # s1 = epg_cpmg(pi, 30, 1000, 200, 10); % IDEAL, CPMG
                echoes[:,z] = np.squeeze(s) # TODO squeeze maybe not necessary


            #% ------------ Slice thickness - -------------------------

            #% ... testes ...
            #% xvec_shortest = 0.8656; (slr_profile.m)
            #% figure, plot(linspace(-xvec_shortest, xvec_shortest, 128), abs(echoes(1,:)))
            #% figure, plot(abs(echoes(1,:))), hold on, plot(abs(echoes(2,:)), 'r')
            #% hold on, plot(abs(echoes(5,:)), 'k')

            # -- Somatorio ao longo da fatia, soma das 128 coordenadas z
            soma[:, i, j]=np.sum(echoes, axis=1)
            # figure; plot(abs(soma(:,:)))


    # =============================================================================
    # %% --- 3 - Variable Dict organized
    # =============================================================================
    D1 = soma
    sz2 = int(D1.shape[1])
    sz3 = int(D1.shape[2])
    Dict = np.zeros((ETL, sz3 * sz2),dtype=float)

    for j in range(0, Dict.shape[1], sz2):
        if sz3 == 1:
            Dict = D1
        else:
            Dict[:, range(j, j + sz2,1)] = D1[:,:, int( (j + sz2) / sz2)]

    return Dict