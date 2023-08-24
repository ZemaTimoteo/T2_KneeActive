#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 28 11:18:25 2021
Perform GRAPPA for pulseq sequences
@author: tfernandes, August 2022
@email: tiagotimoteo@tecnico.ulisboa.pt

Needs:
 -
Toolboxes:
 -

% This code was addapted from MATLAB's code from  Felix Breuer based on opengrappa by Mark Griswold
%
%
%
%  [recon,sigrecon,ws,ws_img,g] = opengrappa(sig,acs,af);
%
%   This is a teaching version of the GRAPPA reconstruction code.
%
%   IN:         sig                reduced data set         (#coils, Ky./af, Kx)
%               acs                autocalibration lines    (#coils, #acs lines, Kx-acslines)
%               af                 Acceleration factor      (integer)
%
%   OUT:        recon              Reconstructed images     (# coils, Ny, Nx)
%               ws                 Weight sets
%
%   Some things to think about when using this code:
%
%           -The ACS lines used for reconstruction are NOT included in the final reconstructed
%           data sets. If you want to do this, you can do it after reconstruction. Please check
%           that they are lined up with the reconstructed data. I have not checked this.
%
%           -Since the ACS lines are not included, feel free to use a different imaging sequence
%           for acquisition of the ACS lines. We have seen some advantage in doing this when
%           the sequence used for the reduced acquisition is very flow sensitive, for example.
%           This can also be faster in many cases.
%
%           -The 4x5 block size is normally optimal for most applications, but if you need
%           something faster, you should try the 4x3 code which uses a smaller block size.
%
%           -The matrix problems here are normally so overdetermined that numerical conditioning
%           is normally not necessary. You could try something out on your problem and see what
%           you get. I have never seen any real advantage.
%
%           -The cyclic boundary condition will not work as programmed for radial and spiral
%           reconstructions. You will have to reprogram this yourself.
%
%           -Finally, the coil index is the first index simply for historical reasons. I will
%           eventually reprogram this. Use the 'permute' command to reorder your matrices if you
%           have the coil index last.
%
%   Please read the license text at the bottom of this program. By using this program, you
%   implicity agree with the license.
%
%   The main points of the license:
%
%   1) This code is strictly for non-commercial applications. The code is protected by
%      multiple patents.
%   2) This code is strictly for research purposes, and should not be used in any
%      diagnostic setting.
%

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
import time
import statistics
import sigpy
import scipy
import pydicom
import twixtools as tt
import ismrmrdtools
import importlib
import nibabel as nib
import tkinter as tk


from grappa import grappa
from tkinter import *
from tkinter import filedialog
from roipoly import RoiPoly
from scipy.io import savemat
from importlib import reload
from statistics import mean
from ismrmrdtools import sense, coils, transform, show, simulation
from nibabel.testing import data_path


sys.path.append('D:/Tiago/Trabalho/2021_2025_PhD/Projects/qMRI_Joint/Code/pythonCode/Reconstruction/T2_EPG')            # Add functions from mytoolbox
sys.path.append('D:/Tiago/Trabalho/2021_2025_PhD/Projects/qMRI_Joint/Code/pythonCode/Toolboxes/mri-sim-py-master/epg')  # Add toolboxes from https://github.com/utcsilab
from dict_pars_generator_seq import dict_pars_generator_seq
from template_match_test import template_match_test
from template_match import template_match
from epgcpmg import FSE_signal

def opengrappa(sig,acs,af,srcx,srcy):

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # ---- 1 - Parameters ----
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Get the size of both the input data and the ref data
    nc     = sig.shape[0]
    ny_red = sig.shape[1]
    nx     = sig.shape[2]

    nc_acs = acs.shape[0]
    nyacs  = acs.shape[1]
    nxacs  = acs.shape[2]
    nslic  = acs.shape[3]
    nechos = acs.shape[4]
    nxacs  = nxacs * nslic * nechos

    ny     = ny_red *af

    if nc_acs!=nc:
        print('Error! The number of coils has to be the same for both inputs!')

    dx = math.floor(srcx/2)
    dy = round((srcy/2-1)*af)

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # ---- 2 - Calculate weights - ---    # %
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    print('Calculating weights ...')

    # % Prepare source and target matrix
    # % Number of source points -> nc*nsry*nsrcx
    # % Number of target points -> nc * (af-1)
    # % number of kernel repetitions in ACS data: (nyacs-(srcy-1)*af)*(nxacs-(srcx-1))

    src = np.zeros([nc*srcy*srcx,(nyacs-(srcy-1)*af)*(nxacs-(srcx-1))], dtype=np.complex64)
    trg = np.zeros([nc*(af-1),(nyacs-(srcy-1)*af)*(nxacs-(srcx-1))], dtype=np.complex64)

    # %   Simple example at af=2 and kernel size 5x4:
    # %
    # %                   -O-O-O-O-O-     1
    # %                   - - - - - -     2
    # %                   -O-O-O-O-O-     3
    # %                   - - -X- - -     4
    # %                   -O-O-O-O-O-     5
    # %                   - - - - - -     6
    # %                   -O-O-O-O-O-     7
    # %                   - - - - - -     8
    # %
    # %   The circles are the source points, and the X are the target points.
    # %   Note: All source points in all coils are used to fit one target point in each individual coil

    cnt = 0             #% This is a very lazy counter. Could be done much faster.
    acs_reshp = acs.reshape(nc,nyacs,nyacs*nechos)
    acs_reshp_test = np.reshape(acs,(nc,nyacs,nyacs*nechos), order='F')
    for xind in range(dx+1-1, nxacs-dx):
        for yind in range(nyacs-(srcy-1)*af):
            cnt
            #% These are the source points (#coils*srcy*srcx)
            aux_vector_acs_for_src = acs_reshp[ : , np.linspace(yind, yind+(srcy-1)*af, num= yind+(srcy-1)*af - yind).astype(int),:]
            vector_acs_for_src     = aux_vector_acs_for_src[ : , : , np.linspace(xind -dx, xind + dx, num= 1 + xind +dx - (xind - dx)).astype(int)]
            src[:,cnt]             = np.squeeze(  vector_acs_for_src.reshape(nc*srcy*srcx, 1, order='F')  )

            #% these are the taget points (#coils*(af-1))
            aux_vector_acs_for_trg = acs_reshp_test[ : , np.linspace(yind+1+dy, yind+dy+af-1, num= 1 + yind+dy+af-1 - (yind+1+dy) ).astype(int),:]
            vector_acs_for_trg     = np.squeeze( aux_vector_acs_for_trg[:,:,xind] )
            trg[:,cnt]             = np.squeeze( vector_acs_for_trg.reshape(nc*(af-1), 1, order='F') )

            cnt = cnt+1

    inv_src = np.linalg.pinv(src)  #% find weights by fitting the source data to target data
    # inv_src = scipy.linalg.pinv(src)
    #
    # Q, R = qr(src)
    # tmp1 = solve(R.T, src.T)
    # tmp2 = solve(R, tmp1)
    # Apinv = tmp2
    ws = np.dot(trg,inv_src)  #% find weights by fitting the source data to target data

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # ---- 2 - Apply weights - ---    # %
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    print('Applyig weights... ')
    #% prepare matrix for convolution
    sigrecon = np.zeros([nc, ny+2*dy+1,nx+2*dx], dtype = np.complex64)

    # % Write undersampled data into zero-padded matrix
    count_dx = 0
    for test_dx in range(dx,sigrecon.shape[2]-dx-1):
        count_dy = 0
        for test_dy in range(dx, ny + dy - 1, 2):
            sigrecon[:,test_dy,test_dx] = sig[:,count_dy,count_dx]
            count_dy                    = count_dy + 1
        count_dx = count_dx + 1

    # %Apply weights to source points
    for xind in range(dx,nx+dx-1):
        for yind in range(0,ny-1,2):
            aux_sigrecon  = sigrecon[ : , np.linspace(yind, yind + (srcy-1)*af, num=  (yind + (srcy-1)*af) - yind).astype(int),:]
            test_sigrecon = aux_sigrecon[ : , : , np.linspace(xind -dx, xind + dx, num= 1 + xind +dx - (xind - dx)).astype(int)]
            src           = np.squeeze( test_sigrecon.reshape(nc*srcy*srcx, 1, order='F') )

            ws_src                                = ws.dot(src)
            sigrecon[:,yind+dy:yind+dy+af-2,xind] = ws_src.reshape(nc, af-1, order='F')

    # %Crop out the good data.
    aux_sigrecon_vf  = sigrecon[ : , np.linspace(dy, ny+dy-1, num=  ny+dy - dy).astype(int),:]
    sigrecon_vf      = aux_sigrecon_vf[ : , : , np.linspace(dx, nx+dx-1, num= nx+dx-1 - dx).astype(int)]

    # recon=ifftshift(ifftshift(ifft(ifft(ifftshift(ifftshift(sigrecon,2),3),[],2),[],3),2),3);
    print('Finish openGrappa... ')

    return sigrecon_vf