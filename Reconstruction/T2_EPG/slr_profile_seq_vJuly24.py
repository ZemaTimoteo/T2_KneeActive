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
from scipy.integrate import cumtrapz


# from get_rf_pulses_shape_vJuly import get_rf_pulses_shape_vJuly

def slr_profile_seq_vJuly24(FA_exc, B1,FA_refoc, ST_exc, ST_refoc,Tsinc,Tsymm, Npoints, dir_rf):


    # =============================================================================
    #%% --- 1 - Set # parameters
    # =============================================================================

    FA_scale = B1
    gamma    = 42.54  # Giromagnetic constant (MHz/T)
    x_ex     = np.linspace(-5,5,Npoints)
    ST_exc   = round(ST_exc,4)
    ST_refoc = round(ST_refoc,4)

    #% Load rf pulses file from pypulseq
    os.chdir(dir_rf)
    FA_exc_degree = int(FA_exc*90/(math.pi/2))
    nameRFexcsave = "rf_pulses_exc{}".format(int(FA_exc_degree)) + "_symm{}".format(Tsymm) + "_st{}".format(ST_exc*1e3) + "mm"
    nameRFexcsave = nameRFexcsave.replace(".", "_")
    rf_exc = scipy.io.loadmat(nameRFexcsave)

    #% Read Excitation waveform
    aux_rf_exct = rf_exc['rf_ex'][0][0]
    aux_t_ex    = rf_exc['t_ex'][0][0]
    aux_G_ex    = rf_exc['G_ex'][0][0]

    # # % Parameters from G, T of excitation pulse and ref excitation pulse
    G_ex    = aux_G_ex/(gamma*1e3) * 1e-1                           # G/cm - intensity of slice selection -> 0.74017 -> G = 7.40 mT/m
    t_ex    = aux_t_ex*1e3                                          # ms - RF pulse duration
    # rf_exct = aux_rf_exct['signal'][0]/gamma * 1e-3 * FA_scale      # mT

    rf_exct = aux_rf_exct['signal'][0]*gamma
    FA_ex   = math.pi/2*FA_scale
    rf_exct = rf_exct/np.sum(rf_exct)*FA_ex   #% Normalize waveform used for excitation and multiply by FA_ex

    # plt.figure()
    # plt.plot(np.abs(mxy))

    [a_ex, b_ex] = abrm(np.transpose(rf_exct), x_ex)
    mxy          = 2*np.conj(a_ex*b_ex)		    # selective excitation

    flip_angle_z = np.zeros( ( 1, mxy.shape[0] ), dtype=float)
    for jj in range(mxy.shape[0]):
        try:
            flip_angle_z[0, jj] = math.asin(abs(mxy[jj]))  # % *(180/pi);
        except:
            flip_angle_z[0, jj] = math.asin(1)  # due to numerical limitations of python
    bb = flip_angle_z[0]
    #plt.figure(), plt.plot(aa), plt.plot(bb)


    # Takes dimensionless x used by abr, and scales it to cm based on a gradient strength g (G/cm) and pulse duration t (ms)
    xvec_ex = x_ex / (4.257 * G_ex * t_ex)    # MATLAB -  xvec_ex = gt2cm(x_ex,G,t) -  # porquê 4.257 - it is from John Pauly, 1992 - Takes dimensionless x used by abr, and scales it to cm based on a gradient strength g (G/cm) and pulse duration t (ms)


    # =============================================================================
    # %% --- 2 -  Echo waveform
    # =============================================================================

    # % ------------ echo 2 --------------------------
    nameRFrefsave = "rf_pulses_ref" + str(FA_refoc)
    nameRFrefsave = nameRFrefsave.replace(".", "_")
    if nameRFrefsave[-2:] == '_0':
        nameRFrefsave = nameRFrefsave[:-2]
    nameRFrefsave = nameRFrefsave + "_sinc{}".format(Tsinc) + "_st{}".format(ST_refoc*1e3) + "mm"
    nameRFrefsave = nameRFrefsave.replace(".", "_")
    rf_refoc = scipy.io.loadmat(nameRFrefsave)
    aux_rf_refoc = rf_refoc['rf_ref'][0][0]
    aux_t_refoc  = rf_refoc['t_ref'][0][0]
    aux_G_refoc  = rf_refoc['G_ref'][0][0]

    # # % Parameters from G and T of excitation pulse and refocusing pulse
    G_rf   = aux_G_refoc/(gamma*1e3) * 1e-1                         # G/cm - 100 G/cm = 1 T/m = 0.1 G/cm = 1 mT/m -> G_rf = 6.16 mT/m - exemplo: 0.61681 G/cm
    t_rf   = aux_t_refoc*1e3                                        # ms
    echo   = aux_rf_refoc['signal'][0]/gamma * 1e-3 *  FA_scale     # mT
    x_rf   = x_ex


    # Get integral to check FArefoc = Gamma Integral ( B1 * dt) - According to Phillips implementation
    echo_2 = aux_rf_refoc['signal'][0]*gamma
    FA_rf  = FA_refoc*(math.pi/180)*FA_scale #   % refoc angle in rad
    echo   = echo_2/sum(echo_2)*FA_rf
    x_rf   = x_ex

    # Get integral to check FAexc = Gamma Integral ( B1 * dt)
    #echo_2_integral = echo * 1e-3  * B1scale                                                       # T
    #timeVect_refoc  = np.linspace(0, t_rf * 1e-3, num=echo_2_integral.shape[0])                    # in s
    #FArefoc_real    = cumtrapz(echo_2_integral * gama_Hz,timeVect_refoc, initial=0) * 180/math.pi  # in degrees

    # FA_ex        = math.pi/2*FA_scale
    # exc_integral = rf_exct/np.sum(rf_exct)*FA_ex #% waveform used for excitation
    # FAexc_real   = cumtrapz(exc_integral * 2*math.pi * gama_Hz,timeVect_exc, initial=0) * 180/math.pi   # in degrees

    #plt.figure()
    #plt.subplot(121)
    #plt.plot(timeVect_refoc*1e3,echo_2_integral*1e6)
    #plt.xlabel('time [ms]')
    #plt.ylabel('RF [uT]')
    #plt.title('RF FA = %5.2f deg' % (FA_refoc))
    #plt.subplot(122)
    #plt.plot(timeVect_refoc*1e3,FArefoc_real)
    #plt.xlabel('time [ms]')
    #plt.ylabel('FA [degree]')
    #plt.title('Gamma Integral[B1dT]')


    [a_rf, b_rf] = abrm(np.transpose(echo), x_rf)
    mxy_se       = 1j*(b_rf*b_rf)   # spin-echo

    flip_angle_z_rf  = np.zeros((1, mxy_se.shape[0]), dtype=float)
    for jj in range(mxy_se.shape[0]):
        flip_angle_z_rf[0,jj] = 2 * math.asin(  math.sqrt( abs( mxy_se[jj] ) )  )   #% *(180/pi);

    xvec_rf = x_rf / (4.257 * G_rf * t_rf) # porquê 4.257 - it is from John Pauly, 1992 - Takes dimensionless x used by abr, and scales it to cm based on a gradient strength g (G/cm) and pulse duration t (ms)


    # Plots
    # plt.figure()
    # plt.subplot(221)
    # plt.plot((rf_exct))
    # plt.subplot(222)
    # plt.plot((echo))
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
    # ind  = np.where(abs( abs(xvec_longest) - xvec_shortest[-1] ) < 2e-2 )       # Garantee that they are aligned
    ind  = np.where(abs( abs(xvec_longest) - xvec_shortest[-1] ) < 0.5 )       # Garantee that they are aligned
    dif  = np.around(abs(  abs(xvec_longest[ind]) - xvec_shortest[-1]  ),10)
    ind2 = np.where(dif == max(dif))

    xvec_longest_cropped = xvec_longest[ ind[0][ind2[0][0]] : ind[0][ind2[-1][-1]] ]
    flip_angle_z_cropped = flip_angle_longest[ ind[0][ind2[0][0]] : ind[0][ind2[-1][-1]]  , 0] #  % n=98 pts

    # ... 3.3 - Interpolate ...
    aux_interp  = np.interp(xvec_shortest, xvec_longest_cropped, flip_angle_z_cropped) # % n=Npoints pts - spline - Quadratic
    test_interp = interp1d(xvec_shortest, aux_interp, kind='cubic')
    flip_angle_z_cropped_interp = test_interp.y

    #plots
    # fig = plt.figure()
    # plt.plot(flip_angle_z_cropped_interp)
    #
    # fig = plt.figure()
    # plt.plot(np.reshape(flip_angle_z_rf,(Npoints,1), order='C'))
    #
    # fig = plt.figure()
    # plt.plot(np.reshape(flip_angle_z,(Npoints,1), order='C'))

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

    # plt.subplot(223)
    # plt.plot((exc_pulse))
    # plt.subplot(224)
    # plt.plot((np.transpose(refoc_pulse)))

    # else:
    #     # ... 3.4 - Get excitation & refoc pulse ...
    #     exc_pulse   = flip_angle_z  # (rad)
    #     refoc_pulse = flip_angle_z_rf



    #os.chdir('D:/Tiago/Trabalho/2021_2025_PhD/Projects/qMRI_Joint/Code/matlabCode/qMRI_tools/Sequence_Optimization/Data/rf_pulses')

    #parametersMATexc = {'exc_pulse': exc_pulse,'original_exc_pulse':rf_exct}
    #nameRFexcsave = "rf_pulses_INTERP_pulseq_exc{}.mat".format(FA_exc)
    #savemat(nameRFexcsave, parametersMATexc)

    #parametersMATref = {'refoc_pulse': refoc_pulse,'original_refoc_pulse':echo}
    #nameRFrefsave = "rf_pulses_INTERP_pulseq_ref{}.mat".format(FA_refoc)
    #savemat(nameRFrefsave, parametersMATref)

    #print("\n\n ----- rf_pulses saved -----  ")

    dataPointsVect = np.linspace(1, Npoints, num=Npoints)  # in s
    FAexc_real     = cumtrapz(exc_pulse*180/math.pi,dataPointsVect, initial=0)      # u.a
    FArefoc_real  = cumtrapz(refoc_pulse[0]*180/math.pi,dataPointsVect, initial=0)   # u.a

    # plt.figure()
    # plt.subplot(121)
    # plt.plot(exc_pulse*180/math.pi)
    # plt.xlabel('Data Points')
    # plt.ylabel('FA [degree]')
    # plt.title('Slice Profile | Area =' + str(np.round(FAexc_real[-1]/Npoints,2)))
    # plt.subplot(122)
    # plt.plot(refoc_pulse[0]*180/math.pi)
    # plt.xlabel('Data Points')
    # plt.ylabel('FA [degree]')
    # plt.title('Slice Profile | Area =' + str(np.round(FArefoc_real[-1]/Npoints,2)))


    return exc_pulse, refoc_pulse