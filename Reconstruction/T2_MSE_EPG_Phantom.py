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
import twixtools as tt
import pathlib as plib
import nibabel as nib
import tkinter as tk
import sys
import math
import statistics
import sigpy
import pydicom
import ismrmrdtools
import importlib
import torch
import json

os.environ["KMP_DUPLICATE_LIB_OK"]="TRUE"


from tkinter import *
from tkinter import filedialog
from roipoly import RoiPoly
from tabulate import tabulate
from scipy.io import savemat
from scipy.io import loadmat
from importlib import reload
from statistics import mean
from ismrmrdtools import sense, coils, transform, show, simulation
from nibabel.testing import data_path


sys.path.append('GitHub/Reconstruction/T2_EPG')            # Add functions from mytoolbox
sys.path.append('GitHub/Toolboxes/mri-sim-py-master/epg')  # Add toolboxes from https://github.com/utcsilab
from dict_pars_generator_seq_vJuly24 import dict_pars_generator_seq_vJuly24
from template_match import template_match
from template_match_test import template_match_test
from grappa import grappa


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
plotTest      = False
saveResults   = True    # Save - 'True' or 'False'
AllLoad       = False    # 'True' Load kspace, image space and coils | 'Fals' Just Load image (much smaller data)
niftiTest     = False   # Read nifti data - 'True' or 'False'
dataType      = 'dat'   # Matlab files: 'mat' OR DICOM files: 'dicom' OR 'dat'
dataSaveJSON  = False

# 2.3 - Settings
maskTest      = False

DictTest      = False
SLR_prof      = True

templateMatch = True
averageTest   = False       # Get average for T2 values within mask

# Sequence Acceleration
accTest = 'LORAKS'   # 'SENSE' OR 'GRAPPA' OR 'LORAKS' OR 'None'

MATLABrunCode = True         # Test Matlab run code


# 2.3 - Directories
dict_dir = "D:/Tiago/Trabalho/2021_2025_PhD/Projects/qMRI_Joint/Data/Dictionaries"
# dir_Data = "D:/Tiago/Trabalho/2021_2025_PhD/Projects/qMRI_Joint/Code/pythonCode/Data_for_code/rf_pulses"
# dir_Data = "D:/Tiago/Trabalho/2021_2025_PhD/Projects/qMRI_Joint/Code/matlabCode/qMRI_tools/Sequence_Optimization/Data/rf_pulses"
dir_Data = "C:/Users/filia/Documents/T2KneeActive/rf_pulses"

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
root           = tk.Tk()
root.withdraw()
file_path      = filedialog.askopenfilename()
dir_set, fname = os.path.split(file_path)
dir_Results    = dir_set + '/results'

# get test
idx_test     = fname.find("MSE_")
end_idx_test = fname.find("_t2")

str_test = fname[idx_test+4:idx_test+7]
if str_test == 'vFA':
    ConstantFA = False  # Test if constant FA 'True' or 'False'

else:
    ConstantFA = True  # Test if constant FA 'True' or 'False'


os.chdir(dir_set)
mat     = scipy.io.loadmat('sequence_info')

    # ... 3.1.2 - load parameters ...
auxTest     = mat['test']  # test
auxTR       = mat['TR']    # Repetition time in (s)
auxTE       = mat['TE']
auxDT       = mat['DT']
auxT2test   = mat['T2test']
auxNslices  = mat['nslices']
auxSt_exc   = mat['st_exc']
auxSt_refoc = mat['st_refoc']
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
auxTsinc    = mat['testSinc']
auxTsymm    = mat['symmetricVal']

test        = auxTest[0][0]
TR          = auxTR[0][0]
TE          = auxTE[0][0]       # Echo Time in (s)
DT          = auxDT[0][0]       # Delta Time in (s)
T2test      = auxT2test[0][0]   # T2 test (s)
TR          = auxTR[0][0]
nslices     = auxNslices[0][0]
st_exc      = auxSt_exc[0][0]
st_refoc    = auxSt_refoc[0][0]
max_grad    = auxMax_grad[0][0]
max_slew    = auxMax_slew[0][0]
flipAngle   = auxFA[0][0]
RFflipAngle = np.transpose(auxRFfa)
gamma       = 42.576e6
nTE         = auxNEchos[0][0]
Nx          = auxNx[0][0]
Ny          = auxNy[0][0]
FOV         = auxFOV[0][0]
duration    = auxDuration[0][0]
delta_k     = auxDelta_k[0][0]
k_width     = auxK_width[0][0]
Tsinc       = auxTsinc[0][0]
Tsymm       = auxTsymm[0][0]

nc          = 1                                        # Number of Coils for Coil compressing
ESP         = TE*1e3                                   # Echo Spacing Time (ms)

    # ... 3.1.3 - set Angles ...
FA_exc    = math.pi/2                                  # Flip-angle Exciatation - in (rad) - along y
phase_exc = math.pi/2                                  # Phase Exciatation - in (rad) - along y

FA_refoc    = RFflipAngle                              # Flip-angle Refocosing - in (degrees) - along x
phase_refoc = np.exp(np.zeros(nTE)/(180/math.pi)*1j)   # Phase Refocosing -  ou zeros(1,nTE);

    # ... 3.1.3 - get Specific GRAPPA parameters ...
if accTest == 'GRAPPA': # +25 & -25 to centre k-space, and the rest (200-50=150) equaly sampled
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

elif accTest == 'LORAKS':
    auxNy_real         = mat['Ny_real']
    auxR               = mat['R']
    auxfullLin         = mat['fullLin']
    auxiniKfull        = mat['iniKfull']
    auxendKfull        = mat['endKfull']
    auxpercFullKspace  = mat['percFullKspace']
    auxmask_Kspace     = mat['mask_Kspace']

    Ny_real        = auxNy_real[0][0]
    R              = auxR[0][0]                     # Acceleration Factor for GRAPPA
    fullLin        = auxfullLin[0][0]               # Full lines in center of k-space
    iniKfull       = Ny_real - auxendKfull[0][0]    # Full lines in center of k-space
    endKfull       = fullLin + iniKfull             # Full lines in center of k-space
    percFullKspace = auxpercFullKspace[0][0]        # Full lines in center of k-space
    mask_Kspace    = auxmask_Kspace          # Mask of k-space



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

    if accTest == 'SENSE':
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


    elif accTest == 'GRAPPA':
    # ... 4.4.1 - Get Reshape K-data ...
        nimgs = Ny_real * nTE * nslices
        [Nx, ncoils, Ny_real_Nechos_nslice_Nreps] = kdata.shape
        aux_k = kdata.reshape(Nx, ncoils, Ny_real, nslices,int(nimgs / Ny_real / nslices))  # [Nx , ncoils , Ny , nslices , nEchoes]

        # save GRAPPA for test in MATLAB
        os.chdir(dir_set)
        GRAPPA_for_data_test = {'aux_k': aux_k, 'R': R, 'iniKfull': iniKfull, 'endKfull': endKfull, 'fullLin': fullLin}
        saveNameFile = "GRAPPA_for_data_test.mat"
        savemat(saveNameFile, GRAPPA_for_data_test)

        if MATLABrunCode: # quit code of python to run grappa_recon_vF.m and get GRAPPA
            print('\n\n WARNING - NEED to Run GRAPPA data on matlab "MSE_preproc_recon_vf.m" \n\n')

    elif accTest == 'LORAKS':
    # ... 4.4.1 - Get Reshape K-data ...
        nimgs = Ny_real * nTE * nslices
        [Nx, ncoils, Ny_real_Nechos_nslice_Nreps] = kdata.shape
        aux_k = kdata.reshape(Nx, ncoils, Ny_real, nslices,
                              int(nimgs / Ny_real / nslices))  # [Nx , ncoils , Ny , nslices , nEchoes]

        # save LORAKS for test in MATLAB
        os.chdir(dir_set)
        saveNameFile = "LORAKS_for_data_test.mat"
        LORAKS_for_data_test = {'aux_k': aux_k, 'R': R, 'iniKfull': iniKfull,
                                'endKfull': endKfull, 'fullLin': fullLin,
                                'mask_Kspace': mask_Kspace}
        savemat(saveNameFile, LORAKS_for_data_test)


        nslices = 1 # TODO - ATENÇAO!!!! - porque nao reconstrui todos os dados por causa de tempo

        # ... 5.3 - Save Json - to use PyMRITools
        if dataSaveJSON:
            # Read path and load data
            path = plib.Path(dir_set)

            k_space = torch.from_numpy(aux_k)
            # the k space in dimension (read, channels, phase, slices, t / echoes)
            # reshape the k'space
            k_space = torch.permute(k_space, (0, 2, 3, 1, 4))

            sampling_mask = torch.from_numpy(mask_Kspace).to(torch.int)
            # the dimensions of your sampling (phase, t / echoes)
            sampling_mask = sampling_mask.unsqueeze(0)  # (1, phase, t / echoes)
            sampling_mask = sampling_mask.expand(k_space.shape[0], -1, -1)
            sampling_mask_tmp = sampling_mask[:, :, None, None, :]\
                                    .expand(-1, -1, k_space.shape[2],k_space.shape[3],-1)\
                                    .to(torch.bool)
            # Toolbox code expects the k space to be 0 filled
            k_space_new = torch.zeros((k_space.shape[0], k_space.shape[0], *k_space.shape[2:]), dtype=k_space.dtype)
            k_space_new[sampling_mask_tmp] = k_space.flatten()

            affine = torch.zeros((4, 4))
            affine[:3, :3] = torch.eye(3)

            # Jochens code is wanting a .pt file of :
            # the k space in dimension (x,y,z,ch,t / echoes)
            # the mask in dimension (x,y,t / echoes)
            # the affine transformation (x,y,t / echoes)
            dir_Json = os.path.join(dir_set,"jsonFiles")
            os.mkdir(dir_Json)
            path_out = plib.Path(dir_Json)
            pathk = path_out.joinpath("jloraks_prep_kspace").with_suffix(".pt")
            torch.save(k_space_new.to(torch.complex64), pathk)
            paths = path_out.joinpath("jloraks_prep_sampling").with_suffix(".pt")
            torch.save(sampling_mask, paths)
            patha = path_out.joinpath("jloraks_prep_affine").with_suffix(".pt")
            torch.save(affine, patha)

            d = {
                "out_path": path.parent.joinpath("results").as_posix(),
                "in_k_space": pathk.as_posix(),
                "in_affine": patha.as_posix(),
                "in_sampling_mask": paths.as_posix(),
                "coil_compression": 32,
                "process_slice": False,
                "radius": 3,
                "c_rank": 150,
                "c_lambda": 0.0,
                # We use the LORAKS S matrix formalism (Jochen claims its working better) And dont use C (i.e. set c_lambda = 0.0).
                # You options for reconstruction are setting the rank and lambda parameters (rank = 150, lambda=0.05, as starting point)
                "s_rank": 200,
                "s_lambda": 0.05,
                "conv_tol": 1e-3,
                "max_num_iter": 15,
                "batch_size": 6,
                "use_gpu": True,
                "gpu_device": 0,
                "visualization": True
                # I always fit one slice at a time, all coils and possibly all echoes onto GPU RAM.
                # If data doesnt fit your GPU RAM, you can try to coil compress
                # (sometimes creates bright and dark spots if taken too far), you can change the batch size,
                # which effectively minimizes the number of echoes reconstructed jointly,
                # you can try to use C instead of S its spawns smaller matrices,
                # or you can try to disable the GPU and use your CPU RAM, at a cost of computation time
            }
            with open(path_out.joinpath("loraks_config").with_suffix(".json"), "w") as jfile:
                json.dump(d, jfile, indent=2)

        if MATLABrunCode: # quit code of python to run grappa_recon_vF.m and get GRAPPA
            print('\n\n WARNING - NEED to Run LORAKS data on matlab "MSE_preproc_recon_vf.m" \n\n')

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


    # --- 5.1 - If SENSE OR GRAPPA OR LORAKS ---
    # https://github.com/ismrmrd/ismrmrd-python-tools/blob/17fe6fbf1bb645112e2ad0023eb4f42afb36421e/parallel_imaging_demo.py
if accTest == 'SENSE':

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
    # tstart = time.time()
    [ coil_map, coil_combined] = coils.calculate_csm_inati_iter(aux_image_csm, smoothing=5, niter=5, thresh=1e-3, verbose=False)
    # print("Inati coil estimation duration: {}s".format(time.time() - tstart))
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
elif accTest == 'GRAPPA':
    if MATLABrunCode:

        os.chdir(dir_set)
        if AllLoad:
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
            mat = scipy.io.loadmat('GRAPPA_recon_matlab_justImag')
            imgMSE_final = mat['recon_Img_full']  # image recon runned in MATLAB for GRAPPA

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

        # ---- Test for LORAKS ---
elif accTest == 'LORAKS':
    if MATLABrunCode:

        os.chdir(dir_set)
        # ... 3.1.2 - load parameters ...
        if AllLoad:
            mat = scipy.io.loadmat('LORAKS_recon_matlab')
            aux_k = mat['kdata_PI']  # kdata runned in MATLAB for GRAPPA
            imgMSE = mat['img_PI']  # image recon runned in MATLAB for GRAPPA

            # ... 5.2.1 - Reconstruction
            imgMSE_final = np.zeros([Nx, Ny, nslices, nTE], dtype=np.complex64)
            for jj in range(nslices):
                for ii in range(nTE):
                    # aux_image     = ifftnd(aux_k[:,:,:,jj,ii], [0,-1])   # ifft
                    aux_image = rms_comb(imgMSE[:, :, :, jj, ii])  # coil compression
                    imgMSE_final[:, :, jj, ii] = aux_image.transpose(1, 0)
                    imgMSE_final = abs(imgMSE_final)
        else:
            mat = scipy.io.loadmat('LORAKS_recon_matlab_justImag')
            imgMSE_final = mat['recon_Img_full']  # image recon runned in MATLAB for GRAPPA

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
        # os.mkdir(dir_Results)
        os.chdir(dir_Results)
        reconData   = {'image': imgMSE_final}
        savemat("reconData.mat", reconData)
        os.chdir(dir_set)


    # ... 5.3 - Plots
if plotTest:
    pointX = 101
    pointY = 165

    # Image space
    plt.figure(figsize=[12, 8])
    plt.title(['Echo', 4 + 1])
    plt.imshow(imgMSE_final[:, :, 0, 4], cmap='gray')
    plt.axis('off')
    plt.clim(np.min(np.min(imgMSE_final)), np.max(np.max(imgMSE_final)))

    for jj in range(nslices):
        # K-space
        # plt.figure(figsize=[12,8])
        # plt.title('k-space')
        # plt.imshow(abs(aux_k[:,1,:,0,1])**0.2, cmap='gray', origin='lower')
        # plt.axis('off')

        # Image space

        plt.figure(figsize=[12, 8])
        for ii in range(nTE):
            plt.subplot(2,round(nTE/2)+1,ii+1)
            plt.title(['Echo', ii+1])
            plt.imshow(imgMSE_final[:,:,jj,ii], cmap='gray')
            plt.axis('off')
            plt.clim(np.min(np.min(imgMSE_final)), np.max(np.max(imgMSE_final)))

    plt.figure(figsize=[12,8])
    for jj in range(nslices):
        # Signal decay
        plt.title('Signal Decay | test %s' % test)
        plt.plot(abs(imgMSE_final[pointX,pointY,jj,:]),label="slice %s" % (jj+1))
    plt.legend()
    plt.show()


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
# parameters
nvials         = 14
#theoric_T2Vial = [939.4, 594.3, 416.5, 267, 184.9, 140.6, 91.76, 64.84, 45.28, 30.62, 19.76, 15.99, 10.47, 8.15]  # test for 1.5T NIST Phantom Lisboa
theoric_T2Vial = [645.8, 423.6, 286.0, 184.8, 134.1, 94.4, 62.5, 45.0, 31.0, 20.1, 15.4, 10.9, 7.6, 5.4]          # test for 3T NIST Phantom Columbia Uni
T2_mask = np.zeros((n_column,n_column,nvials))
if maskTest:

    # ============= MASK FOR T2 ==================================
    for nv in range(nvials):

        # Show the image
        fig = plt.figure()
        plt.imshow(imgMSE_final[:, :, slice_inf, echo],
                   cmap="gist_gray")  # display image 1st echo random slice with good cartilage display
        plt.colorbar()
        plt.title(["draw ROI (mask) for T2 = "+ str(theoric_T2Vial[nv]) + "ms , " + str(nv+1) + " vial in red color"])
        plt.show(block=False)
        # draw new ROI (mask) in red color
        roiMSE_T2 = RoiPoly(color='r', fig=fig)
        # Show the image with the ROI
        fig = plt.figure()
        plt.imshow(imgMSE_final[:, :, slice_inf, echo], cmap="Greys")
        plt.colorbar()
        roiMSE_T2.display_roi()
        # Get mask
        T2_mask[:,:,nv] = roiMSE_T2.get_mask(current_image=imgMSE_final[:, :, slice_inf, echo])

        # Save image mask
        os.chdir(dir_Results)
        plt.savefig('PhantomMask_T2_'+ str(theoric_T2Vial[nv]) + 'ms_vial_' + str(nv+1) + '_slice_1.png')
        print(' --- Performed mask '+ str(nv+1) + 'vial of ' + str(nvials) + ' --- ')

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
    plt.savefig("Phantom_Noise_mask_slice_1.png")

    # save parameters
    if saveResults:
        parametersMASK = {'T2_mask': T2_mask, 'Noise_mask': Noise_mask}
        savemat("roiMSE_T2Phantom_py.mat", parametersMASK)

else:
    os.chdir(dir_Results)
    mask = scipy.io.loadmat('roiMSE_T2Phantom_py')

    aux_T2_mask = mask['T2_mask']
    T2_mask     = np.rot90(np.flip(aux_T2_mask, axis=1))
    # Noise_mask = mask['Noise_mask']

    if plotTest:
        plotT2_mask = T2_mask[:, :, 9]

        fig = plt.figure(2)
        plt.title("T2 Mask")
        plt.imshow(imgMSE_final[:, :, slice_inf, echo], cmap="Greys")
        plt.imshow(plotT2_mask, alpha=.4, cmap="inferno")
        plt.show()

print('\n\n 6 - Sucessfully finished - inputs of Data (Mask) \n\n')




# =============================================================================
#%% --- 6.5 - Read NIFTI data preprocessed
# =============================================================================
if niftiTest:
    nifti_filename = dir_set + "/imagResults_preproc/003_degibbs.nii.gz"

    notPreProc_imgMSE_final = imgMSE_final
    img_nifti               = nib.load(nifti_filename)
    imgMSE_final_nifti      = np.array(img_nifti.dataobj)
    imgMSE_final_nifti      = np.transpose(imgMSE_final_nifti, [1, 0, 2, 3])

    # Figures
    # # imgMSE from .nii
    # fig = plt.figure()
    # plt.subplot(1, 2, 1)
    # plt.title("imgMSE from .nii")
    # plt.imshow((imgMSE_final_nifti[:, :,0,0]), cmap="gist_gray")
    # plt.show
    #
    # # imgMSE from .dat
    # plt.subplot(1, 2, 2)
    # plt.title("imgMSE from .dat")
    # plt.imshow((imgMSE_final[:, :,0,0]), cmap="gist_gray")
    # plt.show

    # Overwrite imgMSE_final
    imgMSE_final_vDat = imgMSE_final
    imgMSE_final      = imgMSE_final_nifti




# =============================================================================
# =============================================================================
#%% --- 7 - Build dictionary -  T2 estimation
# =============================================================================
# =============================================================================

    # ... 7.1 - parameters for dictionary ...
# T1
T1      = 1000                      # ms
# T2
aux1_T2 = np.linspace(1,300,300)   # ms
# aux2_T2 = np.linspace(302,800,250)   # ms
# T2      = np.append(aux1_T2, aux2_T2)
T2      = aux1_T2

# B1 & slr profile
Npoints = 30      # Number of points for The SLR profile (lower as possible)
# B1      = np.linspace(0.8,1.2,21)
# B1      = np.linspace(0.95,1.35,17)
# B1      = np.linspace(0.95,1.25,31)
B1      = np.linspace(0.95,1.50,31)

# Npoints = 23  # Number of points for The SLR profile (lower as possible)
# B1      = np.linspace(0.95,1.35,41)

if DictTest: # Only for one Slice TODO - add more slices

    if SLR_prof:
        FA_refoc = RFflipAngle

    # t0 = time.clock()

    # B1 = np.linspace(0.3,0.7,41)  # test for different B1.
    # B1 = np.linspace(0.93,1.14,22)  # test from B1map

    # T2 = np.linspace(1, 150, 10)  # ms
    # B1 = np.linspace(0.6, 1.4, 2)

        # ... 7.2 - Get dictionary ...  NOTE: 'FA_exc' in 'rad' & 'FA_refoc' in 'degrees'
    [dict_Phantom, pars] = dict_pars_generator_seq_vJuly24(T1=T1,T2=T2,B1=B1,ESP=ESP,ETL=nTE,
                                                   phase_exc=np.array([phase_exc]),
                                                   phase_refoc=phase_refoc,
                                                   FA_exc=np.array([FA_exc]),
                                                   FA_refoc=FA_refoc,
                                                   ST_exc=st_exc,
                                                   ST_refoc=st_refoc,
                                                   Tsinc=Tsinc,
                                                   Tsymm=Tsymm,
                                                   Npoints=Npoints,
                                                   dir_rf = dir_Data,
                                                   SLR=SLR_prof) # 10min

    col_T2 = pars[:, 0]  # indexes of T2
    col_B1 = pars[:, 2] * 100  # indexes of B1

        # ... 7.3 - Normalize dictionary ...
    Dict_Phantom_norm = np.zeros((nTE,dict_Phantom.shape[1]))

    for ii in range(dict_Phantom.shape[1]):
        Dict_Phantom_norm[:,ii] = dict_Phantom[:,ii]/np.linalg.norm(dict_Phantom[:,ii])

    # t1 = time.clock() - t0
    # print("Time elapsed per slice: ", t1)

        # ... 7.4 - Save Dictionary ...
    if saveResults:
        os.chdir(dict_dir)
        parametersDICT = {'Dict_Phantom_norm': Dict_Phantom_norm, 'col_T2': col_T2, 'col_B1': col_B1}
        if ConstantFA:
            if SLR_prof:
                nameMatrix = f"DICT_MSE_CRLB_py_ESP_{round(ESP)}ms_nTE_{nTE}_cfA_T2_{T2test}ms_stEx_{round(st_exc*1e3,2)}mm_stRF_{round(st_refoc*1e3,2)}mm_symm{Tsymm}_sinc{Tsinc}_Npoints{Npoints}_B1min{int(min(B1)*100)}_B1points{B1.shape[0]}_B1max{int(max(B1)*100)}"
            else:
                nameMatrix = f"DICT_dirac_MSE_CRLB_py_ESP_{round(ESP)}ms_nTE_{nTE}_T2_{T2test}ms_stEx_{round(st_exc*1e3,2)}mm_stRF_{round(st_refoc*1e3,2)}mm_symm{Tsymm}_sinc{Tsinc}_Npoints{Npoints}_B1min{int(min(B1)*100)}_B1points{B1.shape[0]}_B1max{int(max(B1)*100)}"
        else:
            if SLR_prof:
                nameMatrix = f"DICT_MSE_CRLB_py_ESP_{round(ESP)}ms_nTE_{nTE}_vfA_T2_{T2test}ms_stEx_{round(st_exc*1e3,2)}mm_stRF_{round(st_refoc*1e3,2)}mm_symm{Tsymm}_sinc{Tsinc}_Npoints{Npoints}_B1min{int(min(B1)*100)}_B1points{B1.shape[0]}_B1max{int(max(B1)*100)}"
            else:
                nameMatrix = f"DICT_dirac_MSE_CRLB_py_ESP_{round(ESP)}ms_nTE_{nTE}_vfA_T2_{T2test}ms_stEx_{round(st_exc*1e3,2)}mm_stRF_{round(st_refoc*1e3,2)}mm_symm{Tsymm}_sinc{Tsinc}_Npoints{Npoints}_B1min{int(min(B1)*100)}_B1points{B1.shape[0]}_B1max{int(max(B1)*100)}"
        nameMatrix = nameMatrix.replace(".", "_")

        savemat(nameMatrix+".mat", parametersDICT)

else:
    os.chdir(dict_dir)
    if ConstantFA:
        if SLR_prof:
            nameMatrix = f"DICT_MSE_CRLB_py_ESP_{round(ESP)}ms_nTE_{nTE}_cfA_T2_{T2test}ms_stEx_{round(st_exc*1e3,2)}mm_stRF_{round(st_refoc*1e3,2)}mm_symm{Tsymm}_sinc{Tsinc}_Npoints{Npoints}_B1min{int(min(B1)*100)}_B1points{B1.shape[0]}_B1max{int(max(B1)*100)}"
            # nameMatrix = f"DICT_B1small_MSE_CRLB_py_ESP_{round(ESP)}ms_nTE_{nTE}_fA_{int(FA_refoc)}"
        else:
            nameMatrix = f"DICT_dirac_MSE_CRLB_py_ESP_{round(ESP)}ms_nTE_{nTE}_cfA_T2_{T2test}ms_stEx_{round(st_exc*1e3,2)}mm_stRF_{round(st_refoc*1e3,2)}mm_symm{Tsymm}_sinc{Tsinc}_Npoints{Npoints}_B1min{int(min(B1)*100)}_B1points{B1.shape[0]}_B1max{int(max(B1)*100)}"
    else:
        if SLR_prof:
            nameMatrix = f"DICT_MSE_CRLB_py_ESP_{round(ESP)}ms_nTE_{nTE}_vfA_T2_{T2test}ms_stEx_{round(st_exc*1e3,2)}mm_stRF_{round(st_refoc*1e3,2)}mm_symm{Tsymm}_sinc{Tsinc}_Npoints{Npoints}_B1min{int(min(B1)*100)}_B1points{B1.shape[0]}_B1max{int(max(B1)*100)}"
        else:
            nameMatrix = f"DICT_dirac_MSE_CRLB_py_ESP_{round(ESP)}ms_nTE_{nTE}_vfA_T2_{T2test}ms_stEx_{round(st_exc*1e3,2)}mm_stRF_{round(st_refoc*1e3,2)}mm_symm{Tsymm}_sinc{Tsinc}_Npoints{Npoints}_B1min{int(min(B1)*100)}_B1points{B1.shape[0]}_B1max{int(max(B1)*100)}"
    nameMatrix = nameMatrix.replace(".", "_")

    Dict = scipy.io.loadmat(nameMatrix)
    Dict_Phantom_norm = Dict['Dict_Phantom_norm'] # TODO alterar retirando o shortTE
    col_T2            = Dict['col_T2']
    col_B1            = Dict['col_B1']

    # if plotTest:
        # fig = plt.figure()
        # plt.plot(abs(Dict_Phantom_norm[:,:]))
        # plt.show


print('\n\n 7 - Sucessfully finished - Build dictionary - T2 estimation \n\n')





# =============================================================================
#%% --- 8 - SNR - calculate Noise
# =============================================================================
if accTest == 'GRAPPA':
    print("With GRAPPA no Noise Mask needed")
elif accTest == 'LORAKS':
    print("With LORAKS no Noise Mask needed")
# else:
#     Noise_mask_reshp = Noise_mask.reshape((1, Nx * Ny), order='C')
#     SNR              = np.zeros((nTE, nvials))
#     SNRlog           = np.zeros((nTE, nvials))
#     aux_Noise_avgX   = np.zeros((nTE, int(np.count_nonzero(Noise_mask_reshp == 1))))
#     Noise_avgX       = np.zeros((nTE))
#
#     X_SNR_MSE = np.squeeze(imgMSE_final)
#     X_SNR_MSE = np.transpose(X_SNR_MSE, (2, 0, 1))
#     X_SNR_MSE = np.reshape(X_SNR_MSE, (nTE, Nx * Ny), order='C')
#     X_SNR     = abs(X_SNR_MSE)
#
#     for ee in range(nTE):
#         nn = 0
#         for i in range(X_SNR.shape[1]):  # SNR obtained at 2nd echo
#             print(i)
#             if Noise_mask_reshp[0,i] == 0:
#                 a=0
#             else:
#                 aux_Noise_avgX[ee, nn] = X_SNR[ee,i]
#                 nn = nn + 1
#
#         Noise_avgX[ee] = np.mean(aux_Noise_avgX[ee, :], axis=0)

# =============================================================================
#%% --- 9 - T2 Matching (w/ Dictionary)
# =============================================================================


# ... 9.1 - initialize parameters ...
nslices          = slice_sup - slice_inf
T2numbValu       = T2_mask.shape[2]

X                = np.zeros((nTE,Nx*Ny,nslices))
T2_dict_map      = np.zeros((Nx,Ny,nslices,T2numbValu))
B1_dict_map      = np.zeros((Nx,Ny,nslices,T2numbValu))
T2_avg_dict      = np.zeros((T2numbValu))
T2_mask_reshp    = np.zeros((T2numbValu,Nx*Ny))
avgX             = np.zeros((nTE,T2numbValu))
StdX             = np.zeros((nTE,T2numbValu))


    # ... 9.2 - Template_matching of the dictionary ...
if templateMatch:
    for aux in range(0, T2numbValu): # cycle in the number of vials
        itt_slice = 0
        # get reshaped mask
        T2_mask_reshp[aux, :] = T2_mask[:, :, aux].reshape((1, Nx * Ny), order='C')

        # if imgMSE_final.shape[2] > 0:  # se mais do que uma slice
        #     imgMSE_final = imgMSE_final[:, :, slice_inf, :]

        if averageTest:
            jj = 0
            aux_avgX = np.zeros((nTE,int(np.count_nonzero(T2_mask_reshp[aux,:]==1))))
        for ind in range(slice_inf,slice_sup,1):
            # t0 = time.perf_counter()

            #X_MSE = np.squeeze(imgMSE_final[:,:,ind,:]) # test for more than 1 slice
            X_MSE = np.squeeze(imgMSE_final[:, :, ind, :])            # test for only 1 slice
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
                        # ind_param[itt_slice,i] = template_match(  Dict_Phantom_norm,  X[:,i,itt_slice] )
                        ind_param[itt_slice,i] = template_match_test(  Dict_Phantom_norm,  X[:,i,itt_slice] , col_T2 , col_B1 )
                        # ind_param[itt_slice,i] = template_match_test(  Dict_Phantom_norm[0:2,:],  X[0:2,i,itt_slice] , col_T2 , col_B1 )

            if averageTest:  # ... 8.2.1 - Get T2 avg dictionary ...
                # get mean & std
                avgX[:,aux] = np.mean(aux_avgX,axis=1)
                StdX[:,aux] = np.std(aux_avgX,axis=1)

                if accTest == 'GRAPPA': # SNR calculate for GRAPPA
                    for ee in range(nTE): # iteration in echoes
                        SNR[ee, aux]    = np.double(avgX[ee, aux] / StdX[ee, aux])
                        SNRlog[ee, aux] = np.double(10 * np.log10(SNR[ee, aux]))

                else: # SNR calculate without GRAPPA
                    for ee in range(nTE):
                        SNR[ee,aux]    = np.double(avgX[ee,aux]/Noise_avgX[ee])
                        SNRlog[ee,aux] = np.double(10*np.log10(SNR[ee,aux]))

                # obtain template match
                ind_param = template_match(Dict_Phantom_norm, avgX[:,aux])
                # ind_param = template_match_test(Dict_Phantom_norm, avgX[aux,:], col_T2, col_B1)

                index_Value      = ind_param
                index_Value      = index_Value.astype(int)
                T2_avg_dict[aux] = col_T2[(0,index_Value-1)]+1

                print("Cartilage T2map - Dictionary, slice: " + str(ind) + ", T2: " + str(theoric_T2Vial[aux]) + ", avg before: " + str(round(T2_avg_dict[aux], 2)) )

            else:  # ... 8.2.2 - Estimate T2 for each point ...
                # t1 = time.perf_counter() - t0
                # print("Time elapsed per slice: ", t1)

                index_Value                     = ind_param[itt_slice,:]
                index_Value                     = index_Value.astype(int)

                if DictTest: # due to the fact that the vector as different dimentions
                    T2_dict = col_T2[(index_Value)]
                    B1_dict = col_B1[(index_Value)]
                else:
                    T2_dict = col_T2[(0,index_Value)]
                    B1_dict = col_B1[(0,index_Value)]

                T2_dict_map [:,:,itt_slice,aux] = np.reshape(T2_dict, (Nx,Ny), order='C') # [Nx, Ny, nSlices,
                B1_dict_map [:,:,itt_slice,aux] = np.reshape(B1_dict, (Nx,Ny), order='C') # [Nx, Ny, nSlices,
                itt_slice = itt_slice + 1


            # ... 9.3. - Save index ...
    if saveResults:
        os.chdir(dir_Results)
        parametersINDEX = {'T2_dict_map': T2_dict_map}
        parametersINDEX_B1 = {'B1_dict_map': B1_dict_map}

        if niftiTest:  # save results for pre_process
            savemat("T2_dictmap_py_preProc.mat", parametersINDEX)
            savemat("B1_dictmap_py_preProc.mat", parametersINDEX_B1)
        else:
            savemat("T2_dictmap_py_allVials.mat", parametersINDEX)
            savemat("B1_dictmap_py_allVials.mat", parametersINDEX_B1)

        if averageTest:
            # ... 9.3.1 - plot avgX ...
            fig = plt.figure()
            plt.plot(avgX)
            plt.show
            # ... 9.3.2 - save T2_avg_dict ...
            T2avg_values_save = {'T2_avg_Vials_dict': T2_avg_dict}
            savemat("T2_dictmap_avg_py.mat", T2avg_values_save)

            # ... 9.3.3 - Table ...
            #prepare variables
            SNR    = np.round( SNR, 2 )
            SNRlog = np.round( SNRlog, 2)
            # StdX   = np.round( StdX, 2)
            percent = np.zeros((nvials))
            for pp in range(nvials):
                percent[pp] = np.round((T2_avg_dict[pp]-theoric_T2Vial[pp])/theoric_T2Vial[pp]*100,2)

            #define header names
            col_names = ["type","vial1", "vial2","vial3", "vial4","vial5", "vial6","vial7", "vial8","vial9", "vial10","vial11", "vial12","vial13", "vial4"]
            #create data
            data = [[   "Theoric", np.double(theoric_T2Vial[0]),theoric_T2Vial[1],theoric_T2Vial[2],theoric_T2Vial[3],theoric_T2Vial[4],
                                        theoric_T2Vial[5],theoric_T2Vial[6],theoric_T2Vial[7],theoric_T2Vial[8],theoric_T2Vial[9],
                                        theoric_T2Vial[10],theoric_T2Vial[11],theoric_T2Vial[12],theoric_T2Vial[13]    ],
                   [   "T2_avg_Dict", T2_avg_dict[0],T2_avg_dict[1],T2_avg_dict[2],T2_avg_dict[3],T2_avg_dict[4],
                                        T2_avg_dict[5],T2_avg_dict[6],T2_avg_dict[7],T2_avg_dict[8],T2_avg_dict[9],
                                        T2_avg_dict[10],T2_avg_dict[11],T2_avg_dict[12],T2_avg_dict[13]     ],
                   # [   "Std_signal", StdX[1,0],StdX[1,1],StdX[1,2],StdX[1,3],StdX[1,4],StdX[1,5],StdX[1,6],StdX[1,7],StdX[1,8],StdX[1,9],
                   #              StdX[1,10],StdX[1,11],StdX[1,12],StdX[1,13]    ],
                   [   "%", percent[0],percent[1],percent[2],percent[3],percent[4],percent[5],percent[6],percent[7],percent[8],percent[9],
                                percent[10],percent[11],percent[12],percent[13]    ],
                   [   "SNR1", SNR[0,0],SNR[0,1],SNR[0,2],SNR[0,3],SNR[0,4],SNR[0,5],SNR[0,6],SNR[0,7],SNR[0,8],SNR[0,9],
                                SNR[0,10],SNR[0,11],SNR[0,12],SNR[0,13]    ],
                   [   "SNR2", SNR[1,0],SNR[1,1],SNR[1,2],SNR[1,3],SNR[1,4],SNR[1,5],SNR[1,6],SNR[1,7],SNR[1,8],SNR[1,9],
                                SNR[1,10],SNR[1,11],SNR[1,12],SNR[1,13]    ],
                   [   "SNR3", SNR[2,0],SNR[2,1],SNR[2,2],SNR[2,3],SNR[2,4],SNR[2,5],SNR[2,6],SNR[2,7],SNR[2,8],SNR[2,9],
                                SNR[2,10],SNR[2,11],SNR[2,12],SNR[2,13]    ],
                   [   "SNR4", SNR[3,0],SNR[3,1],SNR[3,2],SNR[3,3],SNR[3,4],SNR[3,5],SNR[3,6],SNR[3,7],SNR[3,8],SNR[3,9],
                               SNR[3,10],SNR[3,11],SNR[3,12],SNR[3,13]    ],
                   [   "SNR5", SNR[4,0],SNR[4,1],SNR[4,2],SNR[4,3],SNR[4,4],SNR[4,5],SNR[4,6],SNR[4,7],SNR[4,8],SNR[4,9],
                                SNR[4,10],SNR[4,11],SNR[4,12],SNR[4,13]    ],
                   [   "SNR6", SNR[5,0],SNR[5,1],SNR[5,2],SNR[5,3],SNR[5,4],SNR[5,5],SNR[5,6],SNR[5,7],SNR[5,8],SNR[5,9],
                                SNR[5,10],SNR[5,11],SNR[5,12],SNR[5,13]    ],
                   [   "SNR2log", SNRlog[1, 0],SNRlog[1, 1],SNRlog[1, 2],SNRlog[1, 3],SNRlog[1, 4],SNRlog[1, 5],
                                  SNRlog[1, 6],SNRlog[1, 7],SNRlog[1, 8],SNRlog[1, 9],SNRlog[1, 10],SNRlog[1, 11],SNRlog[1,12],SNRlog[1,13]    ]
                    ]
            # display table
            print(tabulate(data, headers=col_names))


else: # load data
    os.chdir(dir_Results)
    if niftiTest:  # save results for pre_process
        Index     = scipy.io.loadmat('T2_dictmap_py_preProc')
        Index     = scipy.io.loadmat('B1_dictmap_py_preProc')
    else:
        Index_T2     = scipy.io.loadmat('T2_dictmap_py_allVials')
        Index_B1     = scipy.io.loadmat('B1_dictmap_py_allVials')
    T2_dict_map = Index_T2['T2_dict_map']
    B1_dict_map = Index_B1['B1_dict_map']

    # ... 9.3 - Plots ...
if plotTest:
    fig = plt.figure()
    plt.imshow(T2_dict_map[:, :, 0,13], cmap="plasma")
    plt.colorbar()
    plt.show

print('\n\n 9 - Sucessfully finished - T2 Dictionary Matching \n\n')





# =============================================================================
#%% --- 10 - Figures - Add T2 map on top of the magnitude T2w images
# =============================================================================

# 10.1 - Figure for T2
itt_slice=0
for ind in range(slice_inf, slice_sup):
    # ... 10.1 - Figures ...
    fig    = plt.figure(ind+62, frameon=False)
    extent = 0, imgMSE_final.shape[0], 0, imgMSE_final.shape[1]

    Z1  = imgMSE_final[:, :, ind, echo]
    Z1  = np.rot90(Z1, k=1, axes=(1, 0))
    im1 = plt.imshow(Z1, cmap=plt.cm.gray, interpolation='nearest',
                     extent=extent)

    if averageTest:
        print("No T2map available")

    else:

        # Get Value and Standard deviation
        avgT2vial = np.zeros(T2numbValu)
        stdT2vial = np.zeros(T2numbValu)

        for aux in range(T2numbValu):
            # ... 10.2 - Mean & Standart Devixation ...
            Z2 = T2_dict_map[:, :, 0, aux]
            Z2[Z2 == 2] = ['nan']

            masked_data = np.ma.masked_array(Z2, np.isnan(Z2))
            # avg_Z2 = np.average(masked_data[np.nonzero(masked_data)])
            # std_Z2 = np.std(masked_data[np.nonzero(masked_data)])
            avg_Z2 = np.average(masked_data)
            std_Z2 = np.std(masked_data)
            T2_Dict_coef = std_Z2 / avg_Z2

            Z2 = np.rot90(Z2, k=1, axes=(1,0))
            im2 = plt.imshow(Z2, cmap=plt.cm.jet, interpolation='nearest', extent=extent,vmin=0, vmax=100)
            avgT2vial[aux] = round(avg_Z2, 2)
            stdT2vial[aux] = round(std_Z2, 2)

            print("Vial T2: " + str(theoric_T2Vial[aux]) + ", avg: " + str(avgT2vial[aux]) + ", std: " + str(
                stdT2vial[aux]))

            # Figures
            #
            #plt.title("Vial T2map - Dictionary, slice: " + str(ind) + ", avg: " + str(round(avg_Z2, 2)) + ", std: " + str(round(std_Z2, 2)))

            # 10.3 - Save
            os.chdir(dir_Results)
            sliNumb = ind+1

    plt.colorbar()
    plt.title('T2 Map (ms)')
    plt.axis('off')
    plt.clim(0, 100)
    plt.show
    if SLR_prof:
        plt.savefig("Cartilage_T2map_Dictionary_slice_%i.png" %sliNumb)
    else:
        plt.savefig("Cartilage_diracSP_T2map_Dictionary_slice_%i.png" %sliNumb)

    itt_slice = itt_slice + 1



# 10.2 - Figure for B1
itt_slice=0
for ind in range(slice_inf, slice_sup):
    # ... 10.1 - Figures ...
    fig3    = plt.figure(ind+81, frameon=False)
    extent1 = 0, imgMSE_final.shape[0], 0, imgMSE_final.shape[1]

    Z11  = imgMSE_final[:, :, ind, echo]
    Z11  = np.rot90(Z11, k=1, axes=(1, 0))
    im11 = plt.imshow(Z11, cmap=plt.cm.gray, interpolation='nearest',
                     extent=extent1)

    if averageTest:
        print("No B1map available")

    else:
        # Get Value and Standard deviation
        avgB1vial = np.zeros(T2numbValu)
        stdB1vial = np.zeros(T2numbValu)

        for aux1 in range(4,T2numbValu):
            # ... 10.2 - Mean & Standart Devixation ...
            Z21 = B1_dict_map[:, :, 0, aux1]
            Z21[Z21 == B1[0]*100] = ['nan']

            masked_data = np.ma.masked_array(Z21, np.isnan(Z21))
            avg_Z21 = np.average(masked_data)
            std_Z21 = np.std(masked_data)
            B1_Dict_coef = std_Z21 / avg_Z21

            Z21 = np.rot90(Z21, k=1, axes=(1,0))
            im21 = plt.imshow(Z21, cmap=plt.cm.jet, interpolation='nearest', extent=extent1, vmin=int(B1[0]*100), vmax=int(B1[-1]*100))

            if not np.ma.is_masked(avg_Z21):
                avgB1vial[aux1] = round(avg_Z21, 2)
                stdB1vial[aux1] = round(std_Z21, 2)
            else:
                avgB1vial[aux1] = np.nan  # or some fallback value
                stdB1vial[aux1] = np.nan


            print("Vial T2: " + str(theoric_T2Vial[aux1]) + ", B1 avg: " + str(avgB1vial[aux1]) + ", B1 std: " + str(
                stdB1vial[aux1]))

            # 10.3 - Save
            os.chdir(dir_Results)
            sliNumb = ind+1

    plt.colorbar()
    plt.title('B1 Map (%)')
    plt.axis('off')
    plt.clim(int(B1[0]*100), int(B1[-1]*100))
    plt.show
    if SLR_prof:
        plt.savefig("Cartilage_B1map_Dictionary_slice_%i.png" %sliNumb)
    else:
        plt.savefig("Cartilage_diracSP_B1map_Dictionary_slice_%i.png" %sliNumb)


    itt_slice = itt_slice + 1

# Save T2 avg & std + B1 avg & std
os.chdir(dir_Results)
allValues = {"theo": theoric_T2Vial, "T2avg": avgT2vial, "T2std": stdT2vial, "B1avg": avgB1vial, "B1std": stdB1vial}
if SLR_prof:
    savemat("avg_&_std_T2_&_B1.mat",allValues)
else:
    savemat("avg_&_std_diracSP_T2_&_B1.mat",allValues)


test
print('\n dir of sequence tested: ', dir_set)
print('\n\n 10 - Sucessfully finished - Add T2 map on top of the magnitude T2w images \n\n')


# =============================================================================
#%% --- 11 - Roi analysis
# =============================================================================


# vectorT2mask = T2_dict_map[idx_phantom_x,idx_phantom_y]
# for ii in range(vectorT2mask.shape[0]):
#     aux_vector += int(vectorT2mask[ii][0])
# mean_roi_dict = mean()
# std_roi_dict  = std2(T2_dict_map(idx_phantom))
#
#     # ... 11.2 - Precision coefficient ...
#