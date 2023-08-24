# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 18:15:25 2021
Definition:
Generation of a dictionary and its indexes (pars) with 3 entries: T2, the
phase of the 90º RF (phi) and the values of B1.
in each entry: the differents states of the magnetizations are present

Functions used:
- slr_profile
- my epg to generate the states over the slice excited (sum) for different T2
- epg.cpmg to generate the states

Inputs:
- T1: constant- relaxation time of the longitudinal magnetization
- T2 : interval- relaxation time of the transversal magnetization
- B1_scale: interval-
- ESP: spacing time between echoes
- ETL: refocusing echoes train length
- refoc_phase = exp(zeros(1,ETL)./(180/pi).*1i); % phase of refoc pulses (0 = about x axis)
- phi = pi/2; % phase of excitation
- FA_exc - Flip angle excitation pulse: vector
- FA_refoc - Flip angle refocusing pulse: vector
- dTE - Spacing Time between Echoes (ms)

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
from slr_profile import slr_profile

def dict_pars_generator(T1,T2,B1,ESP, ETL,refoc_phase,phase_exc,FA_exc,FA_refoc, dir_rf, SLR=True):

    # =============================================================================
    # %% --- 1 - initialization of variable
    # =============================================================================
    refoc_pulse     = []
    exc_pulse       = []
    FA_refocUNIQUE  = np.unique(FA_refoc, 'stable')
    Dict = np.zeros(  (ETL, int((T2.shape)[0] * phase_exc.shape[0] * B1.shape[0]))  )

    # =============================================================================
    # %% --- 2 - Slice profile with SLR method + Dictionary
    # =============================================================================

        # ... 2.1 - Loop over B1 scale values ...
    for jj in range(B1.shape[0]):
        # print(jj)
        # ... 2.2 - slice profile considering SLR profile ...
        if SLR:
            for ii in range(FA_refoc.shape[0]): # only accept one FA_refoc
                [exc_pulse, val_refoc_pulse] = slr_profile(B1[jj], FA_refoc[ii], dir_rf) # corrige com os perfis adaptar para os pulsos da Siemens/pypulseq

            aux_refoc_pulse = val_refoc_pulse[0,:]
            refoc_pulse = np.reshape(  np.repeat(aux_refoc_pulse,ETL), (aux_refoc_pulse.shape[0],ETL) )

            # Dict = Proc_epg(exc_pulse, refoc_pulse, phase_exc, refoc_phase, T1, T2, ESP, ETL)
            # add the entry to LargePhant_dic
            Dict[:, jj * int((T2.shape)[0]) * phase_exc.shape[0]: (jj+1) * int((T2.shape)[0]) * phase_exc.shape[0]] = np.squeeze(  Proc_epg(exc_pulse, refoc_pulse, phase_exc, refoc_phase, T1, T2, ESP, ETL)  ) # add the entry  to LargePhant_dic
            print('\n\n   -> B1: ',jj+1,'/',B1.shape[0], 'test')

        # ... 2.2 - slice profile without considering SLR profile ...
        else:
            exc_pulse   = float(FA_exc)   * B1[jj]
            val_refoc_pulse = float(FA_refoc) * B1[jj]
            refoc_pulse = np.ones(ETL)*val_refoc_pulse

            # Dict = Proc_epg(np.array([exc_pulse]), np.array([refoc_pulse]), np.array([phase_exc]), np.array([refoc_phase]), T1, T2, ESP, ETL)
            Dict[:, jj * int((T2.shape)[0]) * phase_exc.shape[0]: (jj+1) * int((T2.shape)[0]) * phase_exc.shape[0]] = np.squeeze(  Proc_epg(np.array([exc_pulse]), np.array([refoc_pulse]), np.array([phase_exc]), np.array([refoc_phase]), T1, T2, ESP, ETL)   ) # add the entry  to LargePhant_dic


    # =============================================================================
    # %% --- 3 - PARS generation - Variable pars contains the values of T2, the phase of the 90º RF (phi) % and the values of B1.
    # =============================================================================
    #
        # ... 3.1 - initialization    of    variables: ...
    pars = np.zeros( (int((T2.shape)[0]) * phase_exc.shape[0] * B1.shape[0], 3))

        # ... 3.2 - Get variable pars (contains the values of T2, the phase of the 90º RF (phi) and the values of B1) ...
    for j in range(0,pars.shape[0], int((T2.shape)[0]) * phase_exc.shape[0]):
        ind_B1 = int((j +  int((T2.shape)[0]) * phase_exc.shape[0]) / ( int((T2.shape)[0]) * phase_exc.shape[0])) -1

        for i in range(j, j +  int((T2.shape)[0]) * phase_exc.shape[0], int((T2.shape)[0])):
            pars[i: i + int((T2.shape)[0]), 0]  = T2
            ind_phi                                        = int((i + int((T2.shape)[0])) / int((T2.shape)[0]) - ind_B1 *phase_exc.shape[0])-1
            pars[i: i +int((T2.shape)[0]), 1]   = phase_exc[ind_phi] * np.ones(int((T2.shape)[0]))

        pars[j: j + int((T2.shape)[0]) * phase_exc.shape[0], 2] = B1[ind_B1] * np.ones(int((T2.shape)[0]) * phase_exc.shape[0])


    # # =============================================================================
    # # %% --- 4 -  Out parameters
    # # =============================================================================

    return Dict, pars