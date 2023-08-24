# -*- coding: utf-8 -*-
"""
Created on Mon Ago 16 12:24:25 2021
Definition:
Excitation and refocusing Shinnar-LeRoux profiles with slice profile correction (See paper Buonincontri et al., 2015)
in each entry: the differents states of the magnetizations are present

Functions used:
- slr_slice
- my epg to generate the states over the slice excited (sum) for different T2
- epg.cpmg to generate the states

Inputs:
- B1: B1_scale interval
- ESP: spacing time between echoes
- FA_refoc - Flip angle refocusing pulse: vector

@author: tfernandes

Needs:
 -
"""

# =============================================================================
#%% --- 0 - Import functions
# =============================================================================

import os
import scipy
import matplotlib
import tkinter
import numpy as np
matplotlib.use('TKAgg')
import matplotlib.pyplot as plt
import sys
import math
import cmath
import sigpy

from scipy.io import savemat
from tkinter import *
from tkinter import filedialog
from scipy.interpolate import interp1d
from scipy.interpolate import CubicSpline
from sigpy.mri.rf.sim import abrm

def slr_profile(B1,FA_refoc, dir_rf):


    # =============================================================================
    #%% --- 1 - Set # parameters
    # =============================================================================

    FA_scale = B1
    gamma    = 42.54  # Giromagnetic constant (MHz/T)

    # # % Parameters from G and T of excitation pulse and refocusing pulse
    G    = 0.74017    # G/cm - intensity of slice selection -> G = 7.40 mT/m
    t    = 2.9376     # ms - RF pulse duration
    G_rf = 0.61681    # G/cm - 100 G/cm = 1 T/m = 0.1 G/cm = 1 mT/m -> G_rf = 6.16 mT/m
    t_rf = 1.5232     # ms


    # % conversion to cm - perguntar ao Angelo
    # G    = 0.30663  # G/cm - intensity of slice selection
    # t    = 5.3888   # ms - RF pulse duration
    # G_rf = 0.36796  # G/cm
    # t_rf = 5.9840   # ms

    #% Load rf pulses file from pypulseq
    os.chdir(dir_rf)
    rf_pulses = scipy.io.loadmat('rf_pulses')

    #% Read Excitation waveform
    aux_rf_exct = rf_pulses['rf_ex'][0][0]
    rf_exct = aux_rf_exct[0][0]*gamma

    FA_ex   = math.pi/2*FA_scale
    rf_exct = rf_exct/np.sum(rf_exct)*FA_ex #% waveform used for excitation
    x_ex    = np.linspace(-5,5,128)

    [a_ex, b_ex] = abrm(np.transpose(rf_exct), x_ex)
    mxy          = 2*np.conj(a_ex)*b_ex		    # selective excitation

    flip_angle_z = np.zeros( ( 1, mxy.shape[0] ), dtype=float)
    for jj in range(mxy.shape[0]):
        flip_angle_z[0,jj] = math.asin(abs(mxy[jj])) # % *(180/pi);

    # Takes dimensionless x used by abr, and scales it to cm based on a gradient strength g (G/cm) and pulse duration t (ms)
    xvec_ex = x_ex / (4.257 * G * t)    # MATLAB -  xvec_ex = gt2cm(x_ex,G,t)


    # =============================================================================
    # %% --- 2 -  Echo waveform
    # =============================================================================

    # % ------------ echo 2 --------------------------
    aux_rf_refoc = rf_pulses['rf_ref'][0][0]
    echo_2 = aux_rf_refoc[0][0]*gamma

    FA_rf = FA_refoc*(math.pi/180)*FA_scale #   % refoc angle in rad

    #  perguntar ao Angelo
    echo    = echo_2/sum(echo_2)*FA_rf
    x_rf    = x_ex

    [a_rf, b_rf] = abrm(np.transpose(echo), x_rf)
    mxy_se       = 1j*(b_rf*b_rf)   # spin-echo

    flip_angle_z_rf  = np.zeros((1, mxy_se.shape[0]), dtype=float)
    for jj in range(mxy_se.shape[0]):
        flip_angle_z_rf[0,jj] = 2 * math.asin(  math.sqrt( abs( mxy_se[jj] ) )  )   #% *(180/pi);

    xvec_rf = x_rf / (4.257 * G_rf * t_rf)


    # =============================================================================
    # %% --- 3 -  SLICE PROFILE CORRECTION (Buonincontri 2015)
    # =============================================================================

    # ... 3.1 - check the longest profile to crop and interpolate (Exc or Refoc) ...
    if xvec_ex[-1]<xvec_rf[-1]:
         xvec_longest       = xvec_rf
         xvec_shortest      = xvec_ex
         flip_angle_longest = np.transpose(flip_angle_z_rf)
    else:
         xvec_longest       = xvec_ex
         xvec_shortest      = xvec_rf
         flip_angle_longest = np.transpose(flip_angle_z)

    # if flip_angle_longest.shape[0] < xvec_shortest.shape[0]: # interpolate
    # ... 3.2 - check the longest profile to crop and interpolate (Exc or Refoc) ...
    ind  = np.where(abs( abs(xvec_longest) - xvec_shortest[-1] ) < 1e-2 )
    dif  = np.around(abs(  abs(xvec_longest[ind]) - xvec_shortest[-1]  ),10)
    ind2 = np.where(dif == max(dif))

    xvec_longest_cropped = xvec_longest[ ind[0][ind2[0][0]] : ind[0][ind2[-1][-1]] ]
    flip_angle_z_cropped = flip_angle_longest[ ind[0][ind2[0][0]] : ind[0][ind2[-1][-1]]  , 0] #  % n=98 pts

    # ... 3.3 - Interpolate ...
    aux_interp  = np.interp(xvec_shortest, xvec_longest_cropped, flip_angle_z_cropped) # % n=128 pts - spline - Quadratic
    test_interp = interp1d(xvec_shortest, aux_interp, kind='cubic')
    flip_angle_z_cropped_interp = test_interp.y

    #plots
    # fig = plt.figure()
    # plt.plot(flip_angle_z_cropped_interp)
    #
    # fig = plt.figure()
    # plt.plot(np.reshape(flip_angle_z_rf,(128,1), order='C'))
    #
    # fig = plt.figure()
    # plt.plot(np.reshape(flip_angle_z,(128,1), order='C'))

    # ... 3.4 - Get excitation & refoc pulse ...
    if xvec_ex[-1]<xvec_rf[-1]:
        aux_slicThick = np.ndarray.max(flip_angle_z)/2
        aux_Value = abs(aux_slicThick - flip_angle_z[0, 0:int(flip_angle_z.shape[1] / 2)])
        # min()
        exc_pulse   = np.reshape( flip_angle_z ,  ( flip_angle_z.shape[1] )    )
        refoc_pulse = np.reshape( flip_angle_z_cropped_interp , ( 1 ,(flip_angle_z_cropped_interp.shape[0]) )  )   # (rad)
    else:
        exc_pulse   = np.reshape( flip_angle_z_cropped_interp , ( flip_angle_z_cropped_interp.shape[0] )  )    # (rad)
        refoc_pulse = np.reshape(  flip_angle_z_rf  ,  ( 1 , (flip_angle_z_rf.shape[1] )  ) )

    # else:
    #     # ... 3.4 - Get excitation & refoc pulse ...
    #     exc_pulse   = flip_angle_z  # (rad)
    #     refoc_pulse = flip_angle_z_rf

    return exc_pulse, refoc_pulse