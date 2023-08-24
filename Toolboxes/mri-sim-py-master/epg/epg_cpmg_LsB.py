#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Ago 3 19:14:22 2021
Procedure for EPG in CPMG conditions simulations
@author: tfernandes

%	EPG Simulation of CPMG sequence.  First flip angle
%	is 90 about y axis, and others by default are about
%	x-axis (make refoc_pulse complex to change that).
%
%	refoc_pulse = refocusing flip angle or list (radians)
%	etl = echo train length, if single refoc_pulse is to be repeated.
%	T1,T2,esp = relaxation times and echo spacing (arb. units).
%
%	Note that if refoc flip angle is constant, less than pi, and etl > 1 the
%	first refocusing flip angle is the average of 180 and the desired
%	refocusing flip angle, as proposed by Hennig.
%
%	All states are kept, for demo purposes, though this
%	is not really necessary.

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
import cmath

from scipy.io import savemat
from tkinter import *
from tkinter import filedialog
from roipoly import RoiPoly


if PC == 0:
    sys.path.append('D:/Tiago/Trabalho/2021_2025_PhD/Projects/qMRI_Joint/Code/pythonCode/Toolboxes/mri-sim-py-master/epg') # Add toolboxes from https://github.com/utcsilab
    from epgcpmg import relax
    from epgcpmg import grad

def epg_rf_LsB(FpFmZ,alpha,phi):
    """
    Input:
        :param FpFmZ:       3xN vector of F+, F- and Z states.
        :param alpha:   flip angle in radians - abs(refoc_pulse(ech))
        :param phi:     angle of rotation axis from Mx (radians) - angle(refoc_pulse(ech)).

    Output:
        :return: P: or FpFmZ:  Updated FpFmZ state.
                 RR: RF rotation matrix (3x3)
    """

    RR = np.array([
            [ (math.cos(alpha/2))**2, cmath.exp(2*1j*phi)*(math.sin(alpha/2))**2, -1j*cmath.exp(1j*phi)*math.sin(alpha)  ],
            [ cmath.exp(-2*1j*phi)*(math.sin(alpha/2))**2, (math.cos(alpha/2))**2, 1j*cmath.exp(-1j*phi)*math.sin(alpha) ],
            [ -1j/2*cmath.exp(-1j*phi)*math.sin(alpha), 1j/2*cmath.exp(1j*phi)*math.sin(alpha), math.cos(alpha)          ]
         ])
    FpFmZ = np.dot(RR, FpFmZ)

    return FpFmZ, RR

def epg_cpmg_LsB(exc_pulse, exc_phase, refoc_pulse, ETL, T1, T2, dTE, hennig,refoc_phase):
    """
    Input:
        :param exc_pulse:       Magnitude / Profile of excitatory pulse
        :param exc_phase:       Phase of Excitatory Pulse
        :param refoc_pulse:     Magnitude / Profile of Refocusing pulse
        :param ETL:             Echo Train Lenght
        :param T1:              T1 value in s
        :param T2:              T2 value in s
        :param dTE:             Echo Spacing in s
        :param hennig:          1st flip reduced trick (Hennig)
        :param refoc_phase:     Phase of Refocusing Pulse

    Output:
        :return: Matrix Z
    """

    # =============================================================================
    #%% --- 1 - Initial conditions
    # =============================================================================

    FpFmZ      = np.zeros( (3,2*ETL) ) 	    # Allocate all known states, 2 per echo.
    FpFmZ[2,0] = 1                          # Initial condition/equilibrium.

    Pstore = np.zeros( (4*ETL,ETL) )	# Store F,F* states, with 2 gradients per echo
    Zstore = np.zeros( (2*ETL,ETL) )	# Store Z states, with 2 gradients per echo

    # =============================================================================
    #%% --- 2 - 90 excitation
    # =============================================================================
    aux_FpFmZ = epg_rf_LsB(FpFmZ, exc_pulse, exc_phase)  # tip - phase phi in radians
    FpFmZ = aux_FpFmZ[0]
    # =============================================================================
    #%% --- 3 - Refocus pulse
    # =============================================================================

    s = np.zeros((ETL,1))	# Allocate signal vector to store.

    for ech in range( int(ETL) ):
        # -- Left crusher
        FpFmZ = relax(FpFmZ, dTE/2., T1, T2)
        # add diffusion
        FpFmZ = grad(FpFmZ)
        # -- Refoc. RF
        FpFmZ = epg_rf_LsB(FpFmZ, abs(refoc_pulse[ech]), np.angle(refoc_pulse[ech]))
        # -- Right crusher
        FpFmZ = relax(FpFmZ[0], dTE/2., T1, T2) # Need to do FpFmZ[0] because the function epg_rf_LsB gives a tupple
        # add diffusion
        FpFmZ = grad(FpFmZ)

        s[ech] = FpFmZ[(0,0)]                           # Signal is F0 state.

    return s
