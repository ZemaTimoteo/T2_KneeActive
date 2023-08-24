#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 28 11:18:25 2021
Tests T2 MSE estimation with EPG for CRLB T2 protocol optimization in Phantom
@author: tfernandes

Needs:
 - pip install roipoly
 - pip install sigpy
Toolboxes:
 - Add toolbox from https://github.com/ut - mri-sim-py in a toolbox folder
 - Add toolbox from https://github.com/ismrmrd/ismrmrd-python-tools - ISMRMRD Python Tools
 - Add toolbox from https://github.com/pehses/twixtools - Twixtools
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
import time
import statistics
import sigpy
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

if PC == 0:
    sys.path.append('Github/Reconstruction/T2_EPG')            # Add functions from mytoolbox
    sys.path.append('Github/Toolboxes/mri-sim-py-master/epg')  # Add toolboxes from https://github.com/utcsilab
    from dict_pars_generator_seq import dict_pars_generator_seq
    from template_match import template_match
    from epgcpmg import FSE_signal



# =============================================================================
#%% --- 1 - Help Functions
# =============================================================================

def ifftnd(kspace, axes=[-1]):
    from numpy.fft import fftshift, ifftshift, ifftn
    if axes is None:
        axes = range(kspace.ndim)
    img = fftshift(ifftn(ifftshift(kspace, axes=axes), axes=axes), axes=axes)
    img *= np.sqrt(np.prod(np.take(img.shape, axes)))
    return img


def fftnd(img, axes=[-1]):
    from numpy.fft import fftshift, ifftshift, fftn
    if axes is None:
        axes = range(img.ndim)
    kspace = fftshift(fftn(ifftshift(img, axes=axes), axes=axes), axes=axes)
    kspace /= np.sqrt(np.prod(np.take(kspace.shape, axes)))
    return kspace

def rms_comb(sig, axis=1):
    return np.sqrt(np.sum(abs(sig)**2, axis))

print('\n\n 1 - Sucessfully finished - Help Functions \n\n')





# =============================================================================
#%% --- 2 - Settings
# =============================================================================


# 2.1 - Number of slices to study
echo, slice_inf, slice_sup = 1, 1, 2
echo, slice_inf, slice_sup = echo-1, slice_inf-1, slice_sup-1

# 2.2 - Options
plotTest      = True
saveResults   = True     # Save - 'True' or 'Fals'
niftiTest     = False    # Read nifti data - 'True' or 'Fals'
dataType      = 'dat'    # Matlab files: 'mat' OR DICOM files: 'dicom' OR 'dat'

# 2.3 - Settings
maskTest      = False

DictTest      = False
SLR_prof      = False

templateMatch = True
averageTest   = False

SENSEtest     = False        # Test SENSE
GRAPPAtest    = True         # Test GRAPPA
MATLABrunCode = True         # Test Matlab run code

# 2.3 - Directorieshttps://www.evernote.com/shard/s483/nl/86381997/899274d0-f7a5-1d44-bf3a-e5f98a711553?title=A%20Fazer
dict_dir = "Github/Data/Dictionaries"
dir_Data = "Github/Reconstruction/T2_EPG/rf_pulses"

print('\n\n 2 - Sucessfully finished - Settings \n\n')



# =============================================================================
# =============================================================================
#%% --- 3 - Parameters
# =============================================================================
# =============================================================================



# ... part 3.1 - set test ...
type  = 'SCAN'                     # 'PULS' - test102 & test45 or 'SCAN' - test109

# ... part 3.2 - parameters ....
    # ... 3.1.1 - get file ...
root = tk.Tk()
root.withdraw()
file_path = filedialog.askopenfilename()
dir_set, fname = os.path.split(file_path)

if SENSEtest:
    dir_csm   = "C:/Users/filia/Documents/PhD/Projetos/qMRI/Sequences/MSE/test284_GREforSENSE"
    fname_csm = 'meas_MID00690_FID273097_test301v2_MSE_nEchos_noGRAPPA_fov200.dat'
    os.chdir(dir_csm)

os.chdir(dir_set)
mat     = scipy.io.loadmat('sequence_info')

    # ... 3.1.2 - load parameters ...
auxTest     = mat['test']  # test
auxTR       = mat['TR']    # Repetition time in (s)
auxTE       = mat['TE']
auxDT       = mat['DT']
auxNslices  = mat['nslices']
auxSt       = mat['st']
auxMax_grad = mat['max_grad']
auxMax_slew = mat['max_slew']
auxFA       = mat['flipAngle']
auxRFfa     = mat['rf_flipAngle']
auxNEchos   = mat['nEcohs']
auxNx       = mat['Nx']
auxNy       = mat['Ny']
auxFOV      = mat['FOV']
auxDuration = mat['duration']
auxDelta_k  = mat['delta_k']
auxK_width  = mat['k_width']

test        = auxTest[0][0]
TR          = auxTR[0][0]
TE          = auxTE[0][0]       # Echo Time in (s)
DT          = auxDT[0][0]       # Delta Time in (s)
TR          = auxTR[0][0]
nslices     = auxNslices[0][0]
st          = auxSt[0][0]
max_grad    = auxMax_grad[0][0]
max_slew    = auxMax_slew[0][0]
flipAngle   = auxFA[0][0]
RFflipAngle = auxRFfa[0][0]
gamma       = 42.576e6
nTE         = auxNEchos[0][0]
Nx          = auxNx[0][0]
Ny          = auxNy[0][0]
FOV         = auxFOV[0][0]
duration    = auxDuration[0][0]
delta_k     = auxDelta_k[0][0]
k_width     = auxK_width[0][0]

nc          = 1                                        # Number of Coils for Coil compressing
ESP         = TE*1e3                                   # Echo Spacing Time (ms)

    # ... 3.1.3 - set Angles ...
FA_exc    = math.pi/2                                  # Flip-angle Exciatation - in (rad) - along y
phase_exc = math.pi/2                                  # Phase Exciatation - in (rad) - along y

FA_refoc    = RFflipAngle                              # Flip-angle Refocosing - in (degrees) - along x
phase_refoc = np.exp(np.zeros(nTE)/(180/math.pi)*1j)   # Phase Refocosing -  ou zeros(1,nTE);

    # ... 3.1.3 - get Specific SENSE parameters ...
if SENSEtest:
    auxNy_real = mat['Ny_real']
    auxR       = mat['R']

    Ny_real = auxNy_real[0][0]
    R       = auxR[0][0]        # Acceleration Factor for SENSE

    # ... 3.1.4 - get Specific GRAPPA parameters ...
if GRAPPAtest:
    auxNy_real  = mat['Ny_real']
    auxR        = mat['R']
    auxfullLin  = mat['fullLin']
    auxiniKfull = mat['iniKfull']
    auxendKfull = mat['endKfull']

    Ny_real  = auxNy_real[0][0]
    R        = auxR[0][0]           # Acceleration Factor for GRAPPA
    fullLin  = auxfullLin[0][0]     # Full lines in center of k-space
    iniKfull = Ny_real - auxendKfull[0][0]    # Full lines in center of k-space
    endKfull = fullLin + iniKfull  # Full lines in center of k-space

print('\n\n 3 - Sucessfully finished - Define Parametres \n\n')




# =============================================================================
# =============================================================================
#%% --- 4 - Read Data MSE
# =============================================================================
# =============================================================================



if dataType == 'dat':
    # ... 4.1 - get path ...
    # filepath = input("Enter the path of your dicom files: ") # Automatic


    # ... 4.2 - Read .dat file ---
    twix                = tt.read_twix(os.path.join(dir_set, fname))
    image_mdbs          = [mdb for mdb in twix[-1]['mdb'] if mdb.is_image_scan()]   # sort all 'imaging' mdbs into a k-space array ...
    n_line              = 1 + max([mdb.cLin for mdb in image_mdbs])
    n_channel, n_column = image_mdbs[0].data.shape                         # assume that all data were acquired with same number of channels & columns:

    # ... 4.3 - Get K-data ...
    aux_k_data = np.zeros([n_line, n_channel, n_column], dtype=np.complex64)
    for mdb in image_mdbs:
        aux_k_data[mdb.cLin] = mdb.data

    # ... 4.4 - Reshape K-data ...
    # ... 4.4.1 - Get Reshape K-data ...
    kdata   = aux_k_data.transpose(2, 1, 0)

    if SENSEtest:
        nimgs = Ny_real * nTE * nslices
        [Nx,ncoils,Ny_real_Nechos_nslice_Nreps] = kdata.shape
        aux_k   = kdata.reshape(Nx,ncoils,Ny_real,nslices,int(nimgs/Ny_real/nslices))  # [Nx , ncoils , Ny , nslices , nEchoes]

        # auxNx, auxnCoil, auxNy, auxNsli, auxNech = aux_k.shape
        # ... 4.4.2 - get path for GRE CSM...
        # filepath = input("Enter the path of your dicom files: ") # Automatic
        os.chdir(dir_set)

        # ... 4.4.3 - Read .dat file for GRE CSM ...
        twix_csm = tt.read_twix(os.path.join(dir_csm, fname_csm))
        image_mdbs_csm = [mdb for mdb in twix_csm[-1]['mdb'] if mdb.is_image_scan()]  # sort all 'imaging' mdbs into a k-space array ...
        n_line_csm = 1 + max([mdb.cLin for mdb in image_mdbs_csm])

        n_channel_csm, n_column_csm = image_mdbs_csm[0].data.shape  # assume that all data were acquired with same number of channels & columns:

        # ... 4.4.4- Get K-data for GRE CSM ...
        aux_k_data_csm = np.zeros([n_line_csm, n_channel_csm, n_column_csm], dtype=np.complex64)
        for mdb in image_mdbs_csm:
            aux_k_data_csm[mdb.cLin] = mdb.data

        # ... 4.4.5 - Reshape K-data & Get Coil Sensitivity Map (CSM) ...
        kdata_csm = aux_k_data_csm.transpose(2, 1, 0)

        # for GRE with 2 TE
        kdata_csm = kdata_csm[:,:,0::2]
        [Nx_csm, ncoils_csm, Ny_Nechos_nslice_Nreps_csm] = kdata_csm.shape
        Nechos_csm, Nreps_csm, nSlice_csm = 1, 1, 1
        Ny_csm = int(Ny_Nechos_nslice_Nreps_csm/ Nechos_csm / Nreps_csm / nSlice_csm)
        nimgs  = Ny_csm
        aux_k_csm = kdata_csm.reshape(Ny_csm, ncoils_csm, Nx_csm, 1, 1)  # [Nx , ncoils , Ny , nslices , nEchoes]

    elif GRAPPAtest:
        # ... 4.4.1 - Get Reshape K-data ...
        nimgs = Ny_real * nTE * nslices
        [Nx, ncoils, Ny_real_Nechos_nslice_Nreps] = kdata.shape
        aux_k = kdata.reshape(Nx, ncoils, Ny_real, nslices,int(nimgs / Ny_real / nslices))  # [Nx , ncoils , Ny , nslices , nEchoes]

        # save GRAPPA for test in MATLAB
        os.chdir(dir_set)
        GRAPPA_for_data_test = {'aux_k': aux_k, 'R': R, 'iniKfull': iniKfull, 'endKfull': endKfull, 'fullLin': fullLin}
        savemat("GRAPPA_for_data_test.mat", GRAPPA_for_data_test)

        if MATLABrunCode: # quit code of python to run grappa_recon_vF.m and get GRAPPA
            print('\n\n WARNING - NEED to Run GRAPPA data on matlab "grappa_recon_vF.m" \n\n')

    else: # ... 4.4 - Reshape K-data - Normal Acquisition & GRAPPA ...
        # ... 4.4.1 - Get Reshape K-data ...
        nimgs   = Ny * nTE * nslices
        [Nx,ncoils,Ny_Nechos_nslice_Nreps] = kdata.shape
        aux_k   = kdata.reshape(Nx,ncoils,Ny,nslices,int(nimgs/Ny/nslices))  # [Nx , ncoils , Ny , nslices , nEchoes]


print('\n\n 4 - Sucessfully finished - Read Data MSE \n\n')



# =============================================================================
# =============================================================================
#%% --- 5 - load Data MSE & Reconstruction
# =============================================================================
# =============================================================================


    # --- 5.1 - If SENSE OR GRAPPA ---
    # https://github.com/ismrmrd/ismrmrd-python-tools/blob/17fe6fbf1bb645112e2ad0023eb4f42afb36421e/parallel_imaging_demo.py
if SENSEtest:

    ## ======================= Example of Generation of CSM ====================
    # matrix_size = 256
    # csm = simulation.generate_birdcage_sensitivities(matrix_size)
    # phan = simulation.phantom(matrix_size)
    # coil_images = phan[np.newaxis, :, :] * csm
    # show.imshow(abs(coil_images), tile_shape=(4, 2))
    #
    # tstart = time.time()
    # (csm_est2, rho2) = coils.calculate_csm_inati_iter(coil_images)
    # print("Inati coil estimation duration: {}s".format(time.time() - tstart))
    # combined_image2 = np.sum(csm_est2 * coil_images, axis=0)
    #
    # show.imshow(abs(csm_est2), tile_shape=(4, 2), scale=(0, 1))
    # show.imshow(abs(combined_image2), scale=(0, 1))

    ## ======================= Get CSM ====================
    nsli, nimg = 1, 1
    ii, jj = nsli-1, nimg-1
    aux_image_csm = ifftnd(aux_k_csm[:, :, :, 0, 0], [0, -1])  # ifft (Ny,ncoils,Nx)
    aux_image_csm = aux_image_csm.transpose(1,2,0) # (ncoils,Nx,Ny)


    # 5.1.1 - Estimate coil sensitivites from the GRE data
    #calculate_csm_inati_iter(im, smoothing=5, niter=5, thresh=1e-3,verbose=False):
    # Inputs
    #   - im: Input images, [coil, y, x] or [coil, z, y, x].
    #   - smoothing : int or ndarray-like -  Smoothing block size(s) for the spatial axes.
    #   - niter : int  -  Maximal number of iterations to run.
    tstart = time.time()
    [ coil_map, coil_combined] = coils.calculate_csm_inati_iter(aux_image_csm, smoothing=5, niter=5, thresh=1e-3, verbose=False)
    print("Inati coil estimation duration: {}s".format(time.time() - tstart))
    combined_image = np.sum(coil_map * aux_image_csm, axis=0)

    # get coil maps in k-space
    coil_map = coil_map.transpose(2,0,1)   # (ncoils,Nx,Ny) to (Ny,ncoils,Nx)
    k_csm    = fftnd(coil_map, [0, -1])    # ifft (Ny,ncoils,Nx)
    coil_map = coil_map.transpose(1,2,0)   # (ncoils,Nx,Ny)
    k_csm    = k_csm.transpose(1,2,0)      # (ncoils,Nx,Ny)


    if plotTest:
        show.imshow(abs(aux_image_csm), tile_shape=(6, 6))
        show.imshow(abs(coil_map), tile_shape=(6, 6), scale=(0, 1))
        show.imshow(abs(combined_image), scale=(0, 0.001))
        show.imshow(abs(k_csm), tile_shape=(6, 6), scale=(0, 1))

    # ---- Test for SENSE ---
    # 5.1.2 - Reconstruction w/ SENSE
    aux_k = aux_k.transpose(1,0,2,3,4) # transpose data
    for jj in range(nslices):
        # --- Get SENSE map ---
        reload(sense)
        [unmix_sense, gmap_sense] = sense.calculate_sense_unmixing(R, coil_map)
        if plotTest:
            show.imshow(abs(gmap_sense), colorbar=True)

        for ii in range(nTE):
            # --- Reconstruct aliased images ---
            alias_img = transform.transform_kspace_to_image(aux_k[:,:,:,jj,ii], dim=(1, 2)) * np.sqrt(R)
            # --- get Reconstruction w/ SENSE ---
            # recon_sense[:,:,:,jj,ii] = np.squeeze(np.sum(alias_img * unmix_sense, 0))

            if plotTest:
                show.imshow(abs(alias_img))
                show.imshow(abs(recon_sense[:,:,:,jj,ii]), colorbar=True)

    os.chdir(dir_set)
    # save file
    SENSE_for_data_test = {'aux_k_csm': aux_k_csm[:, :, :, 0, 0], 'aux_k': aux_k[:, :, :, jj, ii],
                           'coil_map': coil_map, 'R': R, 'alias_img': alias_img, 'k_csm': k_csm}
    savemat("SENSE_for_data_test.mat", SENSE_for_data_test)

    # ---- Test for GRAPPA ---
if GRAPPAtest:
    if MATLABrunCode:

        os.chdir(dir_set)
        mat = scipy.io.loadmat('GRAPPA_recon_matlab')
        # ... 3.1.2 - load parameters ...
        aux_k  = mat['kdata_PI']  # kdata runned in MATLAB for GRAPPA
        imgMSE = mat['img_PI']    # image recon runned in MATLAB for GRAPPA

        # ... 5.2.1 - Reconstruction
        imgMSE_final = np.zeros([Nx, Ny, nslices,nTE], dtype=np.complex64)
        for jj in range(nslices):
            for ii in range(nTE):
                # aux_image     = ifftnd(aux_k[:,:,:,jj,ii], [0,-1])   # ifft
                aux_image     = rms_comb(imgMSE[:,:,:,jj,ii])                  # coil compression
                imgMSE_final[:,:,jj,ii] = aux_image.transpose(1,0)
                imgMSE_final  = abs(imgMSE_final)


    else:
        af = R
        # K_space for GRAPPA recon
        # aux_k - [Nx , ncoils , Ny , nslices , nEchoes] -
        aux_k2    = aux_k.transpose(2, 0, 1, 3, 4) # aux_k2 - [Nx, Ny, ncoils, nslices, nEchoes]
        orgSizeNx = aux_k2.shape[0]

            # aux for getting points in ky
        ky_1st   = np.linspace(0, iniKfull-1, iniKfull)
        ky_cent  = np.array([iniKfull + af + af*interval for interval in range(round((endKfull-iniKfull)/2))])
        ky_2nd   = np.linspace(endKfull, orgSizeNx-1, orgSizeNx-endKfull)
        ky_grapp = np.concatenate((ky_1st, ky_cent, ky_2nd),axis=0)

        kspac            = np.zeros([aux_k2.shape[1],aux_k2.shape[1],aux_k2.shape[2],aux_k2.shape[3],aux_k2.shape[4]], dtype=np.complex64)
        kspac[::2,:,:,:] = aux_k2[ ky_grapp.astype(int),:,:,:]

        # Center of k-space for calibration
        calib = aux_k2[int(aux_k.shape[0]/2-fullLin/2):int(aux_k.shape[0]/2+fullLin/2),iniKfull:endKfull, :, :].copy()  # [kx, ky, coils]

        # Calibrate a kernel
        kernel_size = (5, 4)

        # reconstruct:
        # ---------------------------
        #     Inputs
        #     -------
        # K-space         : 2D multi-coil k-space data to reconstruct from. Make sure that the missing entries have exact zeros in them.
        # calib           : Calibration data (fully sampled k-space). Coil sensitivity maps (kx, ky, coil).
        # kernel_size     : Size of the 2D GRAPPA kernel (kx, ky).
        # coil_axis       : Dimension holding coil data.  The other two dimensions should be image size: (sx, sy).
        # lamda           : Tikhonov regularization for the kernel calibration.
        # memmap          : Store data in Numpy memmaps.  Use when datasets are too large to store in memory.
        # memmap_filename : Name of memmap to store results in.  File is only saved if memmap=True.
        # silent          : Suppress messages to user.
        #
        #     Returns
        #     -------
        # res             : k-space data where missing entries have been filled in.
        # ---------------------------
        res = grappa(kspac[:,:,:,0,0], calib[:,:,:,0,0], kernel_size, coil_axis=-1, lamda=0.01, memmap=False)

        # Fill out with fully sampled k-space centre ...
        res[int(Nx/2-fullLin/2):int(Nx/2+fullLin/2-1),:,:] = aux_k2[iniKfull+1:endKfull,:,:,0,0]

    #teste 13/05
        aux_reconK= res.transpose(0, 2,1)        # aux_k2 - [Nx, Ny, ncoils, nslices, nEchoes]
    #    aux_image = ifftnd(res, [0, 1])           # ifft
        aux_image = ifftnd(aux_reconK, [0, 1])           # ifft
        aux_image = rms_comb(aux_image)           # coil compression
        imgMSE_final = aux_image.transpose(1, 0)
        imgMSE_final = abs(imgMSE_final)

        plt.figure(figsize=[12,8])
        plt.subplot(2,1,1)
        plt.title('k-space')
        plt.imshow(abs(res[:,:,1])**0.2, cmap='gray', origin='lower')
        plt.axis('off')
        plt.subplot(2,1,2)
        plt.title(['Echo'])
        plt.imshow(imgMSE_final, cmap='gray', origin='lower')
        plt.axis('off')
        plt.clim(np.min(np.min(imgMSE_final)), np.max(np.max(imgMSE_final)))


    # --- 5.2 - Normal Reconstruction ---
else:
    # ... 5.2.1 - Reconstruction
    aux_imgMSE_final = np.zeros([Nx, Ny, nslices,nTE], dtype=np.complex64)
    imgMSE_final = np.zeros([Nx, Ny, nslices,nTE], dtype=np.complex64)
    for jj in range(nslices):
        for ii in range(nTE):
            aux_image     = ifftnd(aux_k[:,:,:,jj,ii], [0,-1])   # ifft
            aux_image     = rms_comb(aux_image)                  # coil compression
            aux_imgMSE_final[:,:,jj,ii] = aux_image.transpose(1,0)
            aux_imgMSE_final            = abs(aux_imgMSE_final)
            imgMSE_final[:,:,jj,ii]     = np.rot90(np.flip(aux_imgMSE_final[:, :, jj, ii], axis=1))
    imgMSE_final  = abs(imgMSE_final)

    # save parameters
    if saveResults:
        dir_Results = dir_set + '/results'
        # os.mkdir(dir_Results)
        os.chdir(dir_Results)
        reconData   = {'image': imgMSE_final}
        savemat("reconData.mat", reconData)
        os.chdir(dir_set)

    # ... 5.3 - Plots
if plotTest:
    plt.figure(figsize=[12,8])
    plt.title('k-space')
    plt.imshow(abs(aux_k[:,1,:,0,1])**0.2, cmap='gray', origin='lower')
    plt.axis('off')

    plt.figure(figsize=[12, 8])
    for ii in range(nTE):
        plt.subplot(2,round(nTE/2)+1,ii+1)
        plt.title(['Echo', ii])
        plt.imshow(imgMSE_final[:,:,slice_inf,ii], cmap='gray')
        plt.axis('off')
        plt.clim(np.min(np.min(imgMSE_final)), np.max(np.max(imgMSE_final)))


print('\n\n 5 - Sucessfully finished - load Data MSE & Reconstruction \n\n')



# =============================================================================
# =============================================================================
#%% --- 6 - inputs of Data (Mask)
# =============================================================================
# =============================================================================



    # ... 6.1 - Create and display images montage ...
# montage_imgMSE = np.array((np.moveaxis(imgMSE_final,-1,2)), dtype=float)     # struct - [Nx, Ny, #echos, #slices]
#
# max_imgMSE = double(max(max(max(abs(montage_imgMSE(:,:, 3,:))))));
# figure;
# montage(montage_imgMSE(:,:, 3,:), [0, max_imgMSE]); % display all slices for 1st echo intensity
#
# montage_imgMSE_disp = permute(imgMSE_final(:,:, slice_inf: slice_sup, 1), [1 2 4 3]); # 1st echo

     # ... 6.2 - segment Phantom T2 ...
if maskTest:
    if nslices == 1:
        # ============= MASK FOR T2 = 8ms ==================================
        # Show the image
        fig = plt.figure()
        plt.imshow(imgMSE_final[:, :, slice_inf, echo],
                   cmap="gist_gray")  # display image 1st echo random slice with good cartilage display
        plt.colorbar()
        plt.title("draw ROI (mask) for T2=8 vial in red color")
        plt.show(block=False)
        # draw new ROI (mask) in red color
        roiMSE_T2_8 = RoiPoly(color='r', fig=fig)
        # Show the image with the ROI
        fig = plt.figure()
        plt.imshow(imgMSE_final[:, :, slice_inf, echo], cmap="Greys")
        plt.colorbar()
        roiMSE_T2_8.display_roi()
        # Get mask
        T2_8_mask = roiMSE_T2_8.get_mask(current_image=imgMSE_final[:, :, slice_inf, echo])
        plt.imshow(T2_8_mask, alpha=.4, cmap="inferno")
        # Save image mask
        os.chdir(dir_set)
        plt.savefig("Phantom_T2_8ms_mask_slice_1.png")

        # ============= MASK FOR T2 = 45ms =================================
        # Show the image
        fig = plt.figure()
        plt.imshow(imgMSE_final[:, :, slice_inf, echo], cmap="gist_gray")  # display image 1st echo random slice with good cartilage display
        plt.colorbar()
        plt.title("draw ROI (mask) for T2=45 vial in red color")
        plt.show(block=False)
        # draw new ROI (mask) in red color
        roiMSE_T2_45 = RoiPoly(color='r', fig=fig)
        # Show the image with the ROI
        fig = plt.figure()
        plt.imshow(imgMSE_final[:, :, slice_inf, echo], cmap="Greys")
        plt.colorbar()
        roiMSE_T2_45.display_roi()
        # Get mask
        T2_45_mask = roiMSE_T2_45.get_mask(current_image=imgMSE_final[:, :, slice_inf, echo])
        plt.imshow(T2_45_mask, alpha=.4, cmap="inferno")
        # Save image mask
        os.chdir(dir_set)
        plt.savefig("Phantom_T2_45ms_mask_slice_1.png")

        # ============= MASK FOR NOISE ==================================
        # Show the image
        fig = plt.figure()
        plt.imshow(imgMSE_final[:, :, slice_inf, echo], cmap="gist_gray")  # display image 1st echo random slice with good cartilage display
        plt.colorbar()
        plt.title("draw ROI (mask) for Noise in red color")
        plt.show(block=False)
        # draw new ROI (mask) in red color
        roiMSE_T2_noise = RoiPoly(color='r', fig=fig)
        # Show the image with the ROI
        fig = plt.figure()
        plt.imshow(imgMSE_final[:, :, slice_inf, echo], cmap="Greys")
        plt.colorbar()
        roiMSE_T2_noise.display_roi()
        # Get mask
        Noise_mask = roiMSE_T2_noise.get_mask(current_image=imgMSE_final[:, :, slice_inf, echo])
        plt.imshow(Noise_mask, alpha=.4, cmap="inferno")
        # Save image mask
        os.chdir(dir_set)
        plt.savefig("Phantom_Noise_mask_slice_1.png")


    else:
        ii = 0
        for ind in range(slice_inf, slice_sup + 1):
            # Show the image
            fig = plt.figure()
            plt.imshow(imgMSE_final[:,:, ind, echo], cmap="gist_gray") # display image 1st echo random slice with good cartilage display
            plt.colorbar()
            plt.title("draw new ROI (mask) in red color")
            plt.show(block=False)
            # draw new ROI (mask) in red color
            roiMSE = RoiPoly(color='r', fig=fig)
            # Show the image with the ROI
            fig = plt.figure()
            plt.imshow(imgMSE_final[:,:, ind, 1], cmap="Greys")
            plt.colorbar()
            roiMSE.display_roi()
            # Get mask
            T2_mask = roiMSE.get_mask(current_image=imgMSE_final[:,:, ind, 1])
            plt.imshow(T2_mask, alpha=.4, cmap="inferno")
            ii = ii + 1
            # Save image mask
            os.chdir(dir_set)
            plt.savefig("CRLB_test_mask_slice_%i.png" % ind)

    os.chdir(dir_set)
    # save parameters
    if saveResults:
        parametersMASK = {'T2_8_mask': T2_8_mask, 'T2_45_mask': T2_45_mask, 'Noise_mask': Noise_mask}
        savemat("roiMSE_T2Phantom_py.mat", parametersMASK)

else:
    # Plot the slice images if plotTest == 'True'
    #if plotTest:

        # if nslices == 1:
        #
        #     fig = plt.figure(1)
        #     plt.imshow(imgMSE_final[:, :, slice_inf, echo], cmap="gist_gray")  # display image 1st echo random slice with good cartilage display
        #     plt.colorbar()
        #     plt.title("T2 Mask")
        #     plt.show(block=False)
        #
        #     os.chdir(dir_set)
        #     mask = scipy.io.loadmat('roiMSE_CRLB_py')
        #     T2_mask = mask['T2_mask']
        #     plt.imshow(T2_mask, alpha=0.3, cmap="inferno")
        #     # Save image mask
        #     os.chdir(dir_set)
        #     plt.savefig("CRLB_test_mask_slice_1.png")
        #
        # else:
        #     for ind in range(slice_inf, slice_sup + 1):
        #
        #         fig = plt.figure(1)
        #         plt.imshow(imgMSE_final[:, :, ind, echo], cmap="gist_gray")  # display image 1st echo random slice with good cartilage display
        #         plt.colorbar()
        #         plt.title("T2 Mask")
        #         plt.show(block=False)
        #
        #         os.chdir(dir_set)
        #         mask = scipy.io.loadmat('roiMSE_CRLB_py')
        #         T2_mask = mask['T2_mask']
        #         plt.imshow(T2_mask, alpha=0.3, cmap="inferno")
        #         # Save image mask
        #         os.chdir(dir_set)
        #         plt.savefig("CRLB_test_mask_slice_%i.png" %ind)

    mask = scipy.io.loadmat('roiMSE_T2Phantom_py')

    # T2_all_mask     = mask['T2_all_mask']

    T2_8_mask     = mask['T2_8_mask']
    T2_45_mask    = mask['T2_45_mask']
    Noise_mask = mask['Noise_mask']

    # T2_mask = np.rot90(np.flip(T2_mask, axis=1))

    if plotTest:
        aux_T2_mask = np.rot90(np.flip(T2_8_mask,axis=1))
        fig = plt.figure(2)
        plt.title("T2 8ms Mask")
        # plt.imshow(T2_mask, cmap="gist_gray")
        plt.imshow(T2_8_mask, cmap="inferno")
        # plt.imshow(np.reshape(T2_mask, (int(np.sqrt(T2_mask.shape[1])),int(np.sqrt(T2_mask.shape[1])))), cmap="Greys")

print('\n\n 6 - Sucessfully finished - inputs of Data (Mask) \n\n')




# =============================================================================
#%% --- 6.5 - Read NIFTI data preprocessed
# =============================================================================
if niftiTest:
    nifti_filename = dir_set + "/imagResults_preproc/003_degibbs.nii.gz"

    notPreProc_imgMSE_final = imgMSE_final
    img                     = nib.load(nifti_filename)
    imgMSE_final             = np.array(img.dataobj)
    imgMSE_final             = np.transpose(imgMSE_final, [1, 0, 2, 3])
    # fig = plt.figure()
    # plt.imshow((imgMSE_final[:, :,0,0]))
    # plt.show

# =============================================================================
#%% --- 7 - Build dictionary - bi-exponencial T2 estimation
# =============================================================================



if DictTest: # Only for one Slice TODO - add more slices

    t0 = time.clock()
        # ... 7.1 - parameters for dictionary ...
    T1 = 600                      # ms
    T2 = np.linspace(1,300,300)   # ms
    B1 = np.linspace(0.6,1.4,81)

    # T2 = np.linspace(1, 150, 10)  # ms
    # B1 = np.linspace(0.6, 1.4, 2)

        # ... 7.2 - Get dictionary ...  NOTE: 'FA_exc' in 'rad' & 'FA_refoc' in 'degrees'
    [dict_Phantom, pars] = dict_pars_generator_seq(T1=T1,T2=T2,B1=B1,ESP=ESP,ETL=nTE,refoc_phase=phase_refoc,phase_exc=np.array([phase_exc]),
                                                   FA_exc=np.array([FA_exc]), FA_refoc=np.array([int(FA_refoc)]),dir_rf = dir_Data, SLR=SLR_prof) # 10min

    col_T2 = pars[:, 0]  # indexes of T2
    col_B1 = pars[:, 2] * 100  # indexes of B1

        # ... 7.3 - Normalize dictionary ...
    Dict_Phantom_norm = np.zeros((nTE,dict_Phantom.shape[1]))

    for ii in range(dict_Phantom.shape[1]):
        Dict_Phantom_norm[:,ii] = dict_Phantom[:,ii]/np.linalg.norm(dict_Phantom[:,ii])

    t1 = time.clock() - t0
    print("Time elapsed per slice: ", t1)

        # ... 7.4 - Save Dictionary ...
    if saveResults:
        os.chdir(dict_dir)
        parametersDICT = {'Dict_Phantom_norm': Dict_Phantom_norm, 'col_T2': col_T2, 'col_B1': col_B1}
        if SLR_prof:
            nameMatrix = f"DICT_MSE_CRLB_py_ESP_{round(ESP)}ms_nTE_{nTE}_fA_{int(FA_refoc)}"
        else:
            nameMatrix = f"DICT_dirac_MSE_CRLB_py_ESP_{round(ESP)}ms_nTE_{nTE}_fA_{int(FA_refoc)}"
        savemat(nameMatrix+".mat", parametersDICT)

else:
    os.chdir(dict_dir)
    if SLR_prof:
        nameMatrix = f"DICT_MSE_CRLB_py_ESP_{round(ESP)}ms_nTE_{nTE}_fA_{int(FA_refoc)}"
    else:
        nameMatrix = f"DICT_dirac_MSE_CRLB_py_ESP_{round(ESP)}ms_nTE_{nTE}_fA_{int(FA_refoc)}"
    Dict = scipy.io.loadmat(nameMatrix)
    Dict_Phantom_norm = Dict['Dict_Phantom_norm'] # TODO alterar retirando o shortTE
    col_T2            = Dict['col_T2']
    col_B1            = Dict['col_B1']

    if plotTest:
        fig = plt.figure()
        plt.plot(abs(Dict_Phantom_norm[:,:]))
        plt.show


print('\n\n 7 - Sucessfully finished - Build dictionary - bi-exponencial T2 estimation \n\n')







# =============================================================================
#%% --- 8 - T2 Matching (w/ Dictionary)
# =============================================================================



    # ... 8.1 - initialize parameters ...
nslices          = slice_sup - slice_inf
T2_test          = np.array([8, 45])
T2numbValu       = T2_test.shape

X                = np.zeros((nTE,Nx*Ny,nslices))
T2_dict_map      = np.zeros((Nx,Ny,nslices,T2numbValu[0]))
T2_mask_reshp    = np.zeros((T2numbValu[0],Nx*Ny))

# T2_mask_reshp{0} = Noise_mask.reshape((1, Nx*Ny), order='C')
T2_mask_reshp[0,:] = T2_8_mask.reshape((1, Nx*Ny), order='C')
T2_mask_reshp[1,:] = T2_45_mask.reshape((1, Nx*Ny), order='C')

# plt.figure()
# plt.plot(T2_mask_reshp)
# plt.show()

    # ... 8.2 - Template_matching of the dictionary ...
if templateMatch:
    for aux in range(T2numbValu[0]):
        itt_slice = 0
        if averageTest:
            jj = 0
            aux_avgX = np.zeros((nTE,int(np.count_nonzero(T2_mask_reshp[aux,:]==1))))
        for ind in range(slice_inf,slice_sup,1):
            t0 = time.clock()

            # X_MSE = np.squeeze(imgMSE_final[:,:,ind,:])
            X_MSE = np.squeeze(imgMSE_final)
            X_MSE = np.transpose(X_MSE, (2, 0, 1))
            X_MSE = np.reshape(X_MSE, (nTE, Nx*Ny), order='C')

            X[:,:,itt_slice] = abs(X_MSE)
            ind_param        = np.zeros((1,X.shape[1]))

            for i in range(X.shape[1]): #47min
                print(i)
                if T2_mask_reshp[aux,i] == 0:
                    ind_param[itt_slice,i] = 1
                else:
                    if averageTest: # get the average time-series for template matching
                        aux_avgX[:,jj] = X[:, i, itt_slice]
                        jj=jj+1
                    else:
                        ind_param[itt_slice,i] = template_match(  Dict_Phantom_norm,  X[:,i,itt_slice] )
                        #ind_param[itt_slice,i] = template_match_test(  Dict_Phantom_norm,  X[:,i,itt_slice] , col_T2 , col_B1 )

            if averageTest:
                avgX      = np.mean(aux_avgX,axis=1)
                ind_param = template_match(Dict_Phantom_norm, avgX, col_T2, col_B1)

                index_Value                     = ind_param
                index_Value                     = index_Value.astype(int)
                T2_dict                         = col_T2[(0,index_Value-1)]+1

                print("Cartilage T2map - Dictionary, slice: " + str(ind) + ", T2: " + str(T2_test[aux]) + ", avg before: " + str(
                    round(T2_dict, 2)) )
            else:
                    # ... 8.2.1 - Estimate of the parameters ...
                t1 = time.clock() - t0
                print("Time elapsed per slice: ", t1)

                index_Value                     = ind_param[itt_slice,:]
                index_Value                     = index_Value.astype(int)
                T2_dict                         = col_T2[(0,index_Value-1)]+1
                T2_dict_map [:,:,itt_slice,aux] = np.reshape(T2_dict, (Nx,Ny), order='C') # [Nx, Ny, nSlices,
                itt_slice = itt_slice + 1

            # ... 8.3. - Save index ...
    if saveResults:
        os.chdir(dir_set)
        parametersINDEX = {'T2_dict_map': T2_dict_map}
        if niftiTest:  # save results for pre_process
            savemat("T2_dictmap_py_preProc.mat", parametersINDEX)
        else:
            savemat("T2_dictmap_py.mat", parametersINDEX)


else: # load data
    os.chdir(dir_set)
    if niftiTest:  # save results for pre_process
        Index     = scipy.io.loadmat('T2_dictmap_py_preProc')
    else:
        Index     = scipy.io.loadmat('T2_dictmap_py')
    T2_dict_map = Index['T2_dict_map']

    # ... 8.3 - Plots ...
if plotTest:
    fig = plt.figure()
    plt.imshow(T2_dict_map[:, :, 0,1], cmap="plasma")
    plt.colorbar()
    plt.show

print('\n\n 8 - Sucessfully finished - T2 Dictionary Matching \n\n')


# =============================================================================
#%% --- 9 - Figures - Add T2 map on top of the magnitude T2w images
# =============================================================================

itt_slice=0

for ind in range(slice_inf, slice_sup):
    # ... 9.1 - Figures ...
    fig    = plt.figure(ind+61, frameon=False)
    extent = 0, imgMSE_final.shape[0], 0, imgMSE_final.shape[1]

    # Z1  = imgMSE_final[:, :, ind, echo]
    Z1  = imgMSE_final[:, :, ind, echo]
    im1 = plt.imshow(Z1, cmap=plt.cm.gray, interpolation='nearest',
                     extent=extent)

    for aux in range(T2numbValu[0]):
        # ... 9.2 - Mean & Standart Devixation ...
        Z2        = T2_dict_map[:, :, itt_slice,aux]
        Z2[Z2==2] = ['nan']

        masked_data  = np.ma.masked_array(Z2, np.isnan(Z2))
        avg_Z2       = np.average(masked_data)
        std_Z2       = np.std(masked_data)
        T2_Dict_coef = std_Z2 / avg_Z2

        im2 = plt.imshow(Z2, cmap=plt.cm.jet, interpolation='nearest',
                         extent=extent)
        print("Cartilage T2map - Dictionary, slice: " +str(ind) + ", T2: " +str(T2_test[aux]) +", avg: "+str(round(avg_Z2, 2)) + ", std: " + str(round(std_Z2, 2)))

    plt.colorbar()
    plt.title("Cartilage T2map - Dictionary, slice: " +str(ind) + ", avg: "+str(round(avg_Z2, 2)) + ", std: " + str(round(std_Z2, 2)) )
    plt.show
    plt.clim(0,48)

    # 9.3 - Save
    os.chdir(dir_set)
    sliNumb = ind+1
    plt.savefig("Cartilage_T2map_Dictionary_slice_%i.png" %sliNumb)

    itt_slice = itt_slice + 1


print('\n\n 9 - Sucessfully finished - Add T2 map on top of the magnitude T2w images \n\n')


# =============================================================================
#%% --- 10 - Roi analysis
# =============================================================================


# vectorT2mask = T2_dict_map[idx_phantom_x,idx_phantom_y]
# for ii in range(vectorT2mask.shape[0]):
#     aux_vector += int(vectorT2mask[ii][0])
# mean_roi_dict = mean()
# std_roi_dict  = std2(T2_dict_map(idx_phantom))
#
#     # ... 6.2 - Precision coefficient ...
#