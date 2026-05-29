
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
import statistics
import sigpy
import pydicom
import twixtools as tt
import ismrmrdtools
import importlib
import nibabel as nib
import tkinter as tk
import winsound
import pylab
import nibabel as nib
import shutil

os.environ["KMP_DUPLICATE_LIB_OK"]="TRUE"

from pathlib import Path
from tkinter import *
from tkinter import filedialog
from roipoly import RoiPoly
from tabulate import tabulate
from scipy.io import savemat
from importlib import reload
from statistics import mean
from nibabel.testing import data_path
import SimpleITK as sitk

if PC == 1:

    os.chdir('/home/tfernandes/Documents/Projetos/Project_lfStroke/Code/1_pythonCodes/Spiral/pulses')

if PC == 0:
    sys.path.append('Github/Reconstruction/T2_EPG')            # Add functions from mytoolbox
    sys.path.append('Github/Toolboxes/mri-sim-py-master/epg')  # Add toolboxes from https://github.com/utcsilab
    from dict_pars_generator_seq_vJuly24 import dict_pars_generator_seq_vJuly24
    from dict_pars_generator_seq_vJuly25 import dict_pars_generator_seq_vJuly25
    from template_match_test import template_match_test
    from template_match_vComplex import template_match_vComplex
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

def rms_comb(sig, axis=1): # root mean square - coil combination
    return np.sqrt(np.sum(abs(sig)**2, axis))

print('\n\n 1 - Sucessfully finished - Help Functions \n\n')





# =============================================================================
#%% --- 2 - Settings
# =============================================================================


# 2.2 - Options
plotTest      = False
finalPlots    = False
saveResults   = True    # Save - 'True' or 'False'
niftiTest     = True   # Read nifti data - 'True' or 'False'
dataType      = 'noLOAD'   # Matlab files: 'mat' OR DICOM files: 'dicom' OR 'dat' OR no LOAD 'noLOAD'

# 2.3 - Settings
SLR_prof       = True
ResampledTest  = False
testCrop       = False
complexDict    = False
DenoiseTest    = False

# Sequence Acceleration
accTest     = 'LORAKS'         # 'GRAPPA' OR 'LORAKS' OR 'noAcc' - No Acceleration
str_test    = 'vFA'            # Type of test - 'cFA' OR 'vFA'
testSegment = 'AutoSegm'       # Type of segmentation: 'AutoSegm' - automatic Segmentation OR 'ManualSegm' - Manual Segmentation (draw in ITK-Snap and saved in nifti)


MATLABrunCode = True         # Test Matlab run code
maskTest      = True
HistogramTest = True

if ResampledTest:
    strResampl = '_Resampled'
else:
    strResampl = ''


# 2.3 - Directories
dict_dir                 = "Github/Data/Dictionaries"
dir_Data                 = "Github/Data/rf_pulses"
base_segmentation_folder = "Github/Data/Segmented"

print('\n\n 2 - Sucessfully finished - Settings \n\n')



# =============================================================================
# =============================================================================
#%% --- 3 - Parameters
# =============================================================================
# =============================================================================



# ... part 3.1 - set test ...
type  = 'SCAN'                     # 'PULS' - test102 & test45 or 'SCAN' - test109

# ... part 3.2 - parameters ...
    # ... 3.1.1 - get file ...
root           = tk.Tk()
root.withdraw()
file_path  = filedialog.askdirectory(title="Select DataAcquisition Folder")

# cycle through all folders
parent_dir = file_path
subjects = [f"Subject{str(i).zfill(2)}" for i in range(1, 2)]
days_reps = ['Day1_Rep1', 'Day1_Rep2', 'Day2_Rep1', 'Day2_Rep2']
# days_reps = ['Day2_Rep1']
test_folders = ['Test731', 'Test732', 'Test733']

for subject in subjects:
    subject_path = os.path.join(parent_dir, subject)

    for dr in days_reps:
        dr_path = Path(os.path.join(subject_path, dr)).as_posix()

        # 🔍 Dynamically find folders that start with 'Test'
        test_folders = [
            name for name in os.listdir(dr_path)
            if os.path.isdir(os.path.join(dr_path, name)) and name.startswith(("test731", "test732", "test733"))
            # if os.path.isdir(os.path.join(dr_path, name)) and name.startswith(("test733"))
        ]

        for test_folder in test_folders:
            dir_set = Path(os.path.join(dr_path, test_folder)).as_posix()

            # Continue to pass through: Test732 | subj4 | day2 | rep1 - due to misssing data
            if dr == 'Day2_Rep1' and subject == 4 and test_folder.startswith('test732'):
                continue  # Skip this iteration


            # Find the .dat file
            for file in os.listdir(dir_set):
                if file.endswith('.dat'):
                    fname = Path(os.path.join(dir_set, file)).as_posix()

            # Settings
            testSeq     = test_folder[4:7]      # Define Test number in variable
            dir_Results = dir_set + '/results' + strResampl
            maskFld     = dir_set + '/segmentation'
            if not os.path.exists(maskFld):
                os.makedirs(maskFld)

            # get test
            idx_test     = fname.find("MSE_")
            end_idx_test = fname.find("_t2")

            # str_test = fname[idx_test+4:idx_test+7]

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

                # ... 3.1.2 - set Angles ...
            FA_exc    = math.pi/2                                  # Flip-angle Exciatation - in (rad) - along y
            phase_exc = math.pi/2                                  # Phase Exciatation - in (rad) - along y

            FA_refoc    = RFflipAngle                              # Flip-angle Refocosing - in (degrees) - along x
            phase_refoc = np.exp(np.zeros(nTE)/(180/math.pi)*1j)   # Phase Refocosing -  ou zeros(1,nTE);

                # ... 3.1.3 - get Specific parameters ...
            if accTest == 'GRAPPA': # +25 & -25 to centre k-space, and the rest (200-50=150) equaly sampled
                auxNy_real  = mat['Ny_real']
                auxR        = mat['R']
                auxfullLin  = mat['fullLin']
                auxiniKfull = mat['iniKfull']
                auxendKfull = mat['endKfull']

                Ny_real  = auxNy_real[0][0]
                R        = auxR[0][0]                       # Acceleration Factor for GRAPPA
                fullLin  = auxfullLin[0][0]                 # Full lines in center of k-space
                iniKfull = Ny_real - auxendKfull[0][0]      # Full lines in center of k-space
                endKfull = fullLin + iniKfull               # Full lines in center of k-space

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
                mask_Kspace    = auxmask_Kspace                 # Mask of k-space

            elif accTest == 'noAcc':
                auxNy         = mat['Ny']
                Ny            = auxNy[0][0]



            # ... part 3.3 - B1map ....
            B1test = 'noB1pre'

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

                if accTest == 'GRAPPA':
                # ... 4.4.1 - Get Reshape K-data ...
                    nimgs = Ny_real * nTE * nslices
                    [Nx, ncoils, Ny_real_Nechos_nslice_Nreps] = kdata.shape
                    aux_k = kdata.reshape(Nx, ncoils, Ny_real, nslices,int(nimgs / Ny_real / nslices))  # [Nx , ncoils , Ny , nslices , nEchoes]

                    # save GRAPPA for test in MATLAB
                    os.chdir(dir_set)
                    GRAPPA_for_data_test = {'aux_k': aux_k, 'R': R, 'iniKfull': iniKfull, 'endKfull': endKfull, 'fullLin': fullLin}
                    savemat("GRAPPA_for_data_test.mat", GRAPPA_for_data_test)

                    if MATLABrunCode: # quit code of python to run grappa_recon_vF.m and get GRAPPA
                        print('\n\n WARNING - NEED to Run GRAPPA data on matlab "MSE_preproc_recon_vf.m" \n\n')

                    nslices = 1 # TODO - ATENÇAO!!!! - porque nao reconstrui todos os dados por causa de tempo

                elif accTest == 'LORAKS':
                # ... 4.4.1 - Get Reshape K-data ...
                    nimgs = Ny_real * nTE * nslices
                    [Nx, ncoils, Ny_real_Nechos_nslice_Nreps] = kdata.shape
                    aux_k = kdata.reshape(Nx, ncoils, Ny_real, nslices,int(nimgs / Ny_real / nslices))  # [Nx , ncoils , Ny , nslices , nEchoes]

                    # save LORAKS for test in MATLAB
                    os.chdir(dir_set)
                    LORAKS_for_data_test = {'aux_k': aux_k, 'R': R, 'iniKfull': iniKfull, 'endKfull': endKfull, 'fullLin': fullLin, 'mask_Kspace': mask_Kspace}
                    savemat("LORAKS_for_data_test.mat", LORAKS_for_data_test)


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


                # --- 5.1 - If GRAPPA OR LORAKS ---
                # https://github.com/ismrmrd/ismrmrd-python-tools/blob/17fe6fbf1bb645112e2ad0023eb4f42afb36421e/parallel_imaging_demo.py
            if accTest == 'GRAPPA':
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
                            imgMSE_final            = abs(imgMSE_final)


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
                    if complexDict:
                        if DenoiseTest:
                            mat = scipy.io.loadmat('LORAKS_recon_matlab_CC_Denoise')
                        else:
                            mat = scipy.io.loadmat('LORAKS_recon_matlab_CoilComplex')
                        aux_imgMSE_final = mat['recon_Img_fullComplx']  # image recon runned in MATLAB for GRAPPA

                    else:
                        if DenoiseTest:
                            mat = scipy.io.loadmat('LORAKS_recon_matlab_justImag_Denoise')

                        else:
                            mat = scipy.io.loadmat('LORAKS_recon_matlab_justImag')

                        # aux_k = mat['kdata_PI']  # kdata runned in MATLAB for GRAPPA
                        aux_imgMSE_final = mat['recon_Img_full']  # image recon runned in MATLAB for GRAPPA

                    [nx, ny, nsl, netl] = aux_imgMSE_final.shape
                    imgMSE_final = np.zeros([Nx, Ny, nslices, nTE])

                    # Guarantee that the image size is the expected
                    if netl != nTE:
                        imgMSE_final   = aux_imgMSE_final[:, :, :, 0:nTE]
                        recon_Img_full = {'recon_Img_full': imgMSE_final}
                        savemat("LORAKS_recon_matlab_justImag.mat", recon_Img_full)
                    else:
                        imgMSE_final   = aux_imgMSE_final
                        recon_Img_full = []

                    del aux_imgMSE_final, nx, ny, nsl, netl, recon_Img_full  # deletes variables


                    if testSeq == '731' or testSeq == '732':
                        aux_imgMSE_final = np.rot90(imgMSE_final[:, ::-1, :, :])
                    else:
                        aux_imgMSE_final = np.rot90(imgMSE_final[:, ::-1, :, :])

                        if subject == 'Subject01' or subject == 'Subject02':
                            if  testSeq == '733' and dr == 'Day1_Rep1':
                                aux_imgMSE_final = np.rot90(imgMSE_final[:, ::-1, :, :])
                                aux_imgMSE_final = np.rot90 (aux_imgMSE_final[:, ::-1, :, :])

                    imgMSE_final = aux_imgMSE_final

                # --- 5.2 - Normal Reconstruction ---
            else:
                # ... 5.2.1 - Reconstruction
                aux_imgMSE_final = np.zeros([Nx, Ny, nslices,nTE], dtype=np.complex64)
                imgMSE_final = np.zeros([Nx, Ny, nslices,nTE], dtype=np.complex64)

                if saveResults:
                    if os.path.exists(dir_Results):
                        os.chdir(dir_Results)
                    else:
                        os.mkdir(dir_Results)
                        os.chdir(dir_Results)

                for jj in range(nslices):
                    for ii in range(nTE):
                        aux_image     = ifftnd(aux_k[:,:,:,jj,ii], [0,-1])   # ifft
                        aux_image     = rms_comb(aux_image)                  # coil compression
                        aux_imgMSE_final[:,:,jj,ii] = aux_image.transpose(1,0)
                        aux_imgMSE_final            = abs(aux_imgMSE_final)
                        imgMSE_final[:,:,jj,ii]     = np.rot90(aux_imgMSE_final[:, :, jj, ii])

                        # Save as nifti
                        # Create a NIfTI image object
                    if saveResults:
                        os.chdir(dir_Results)
                        nifti_image = nib.Nifti1Image(imgMSE_final[:,:,jj,:], affine=np.eye(4))  # Identity matrix as affine
                        # Save the image to a NIfTI file
                        nameFile = 'Recon_' + str(jj+1) +'slice_AllEcho.nii'
                        nib.save(nifti_image, nameFile)
                        os.chdir(dir_set)
                imgMSE_final  = abs(imgMSE_final)

                # save parameters
                if saveResults:
                    reconData   = {'image': imgMSE_final}
                    KspaceData  = {'k_space': aux_k}
                    savemat("reconData.mat", reconData)
                    savemat("kspaceData.mat", KspaceData)
                    os.chdir(dir_set)




                # ... 5.3 - Plots
            if plotTest:
                plt.figure(figsize=[12,8])
                plt.title('k-space')
                plt.imshow(abs(aux_k[:,1,:,0,1])**0.2, cmap='gray', origin='lower')
                plt.axis('off')

                # Plot per echo
                plt.figure(figsize=[12, 8])
                for ii in range(nTE):
                    plt.subplot(2,round(nTE/2)+1,ii+1)
                    plt.title(['Echo', ii])
                    plt.imshow(np.abs(imgMSE_final[:,:,5,ii]), cmap='gray')
                    plt.axis('off')
                    plt.clim(np.min(np.min(np.abs(imgMSE_final))), np.max(np.max(np.abs(imgMSE_final))))

                # Plot per slice
                plt.figure(figsize=[12, 8])
                for jj in range(nslices):
                    plt.subplot(2,round(nslices/2)+1,jj+1)
                    plt.title(['slice', jj])
                    plt.imshow(np.abs(imgMSE_final[:,:,jj,0]), cmap='gray')
                    plt.axis('off')
                    plt.clim(np.min(np.min(np.abs(imgMSE_final))), np.max(np.max(np.abs(imgMSE_final))))

            print('\n\n 5 - Sucessfully finished - load Data MSE & Reconstruction \n\n')



            # =============================================================================
            # =============================================================================
            #%% --- 6 - inputs of Data (Mask)
            # =============================================================================
            # =============================================================================


                 # ... 6.1 - segment Knee - Automatic Cartilage Segmentation T2 ...

            if testSegment == 'AutoSegm':


                # Build the expected folder name
                # folder_name_pattern = f"{subject}_{dr}_{testSeq}_Right"
                folder_name_pattern = f"{subject}_{dr}_{testSeq}.nii.gz"

                # List all folders in the segmentation directory
                all_folders = os.listdir(base_segmentation_folder)

                # Find matching folder
                matched_folder = next((f for f in all_folders if f == folder_name_pattern), None)

                if matched_folder:
                    # src_file = os.path.join(base_segmentation_folder, matched_folder, "Split_labels_Native.nii.gz")
                    src_file = os.path.join(base_segmentation_folder, matched_folder)
                    dst_file = os.path.join(maskFld, "Segmentation.nii.gz")

                    if os.path.exists(src_file):
                        shutil.copyfile(src_file, dst_file)
                        print(f"Copied to: {dst_file}")
                    else:
                        print("Source file does not exist:", src_file)
                else:
                    print("No matching folder found for:", folder_name_pattern)


                os.chdir(maskFld)
                filename = 'Segmentation.nii.gz'
                aux_knee_mask = nib.load(filename)
                T2_mask = aux_knee_mask.get_fdata()
                aux_T2_mask = aux_knee_mask.get_fdata()

                # Label of segmentation
                labels = np.arange(44)+1 # according to Amir, 2022
                T2_mask        = np.zeros(aux_T2_mask.shape)
                T2_mask_full   = np.zeros(aux_T2_mask.shape)
                T2_mask_single = np.zeros(aux_T2_mask.shape)

                for jj in range(nslices):
                    aux2_T2_mask_full   = np.zeros((aux_T2_mask.shape[0],aux_T2_mask.shape[1]))
                    aux2_T2_mask_single = np.zeros((aux_T2_mask.shape[0],aux_T2_mask.shape[1]))
                    aux_mask     = aux_T2_mask[:, :, jj]
                    mask         = np.isin(aux_mask, labels)

                    # each with each condition
                    aux2_T2_mask_full[mask]   = aux_mask[mask]
                    # all to 1
                    aux2_T2_mask_single[mask] = 1


                    # Test to use - use single
                    if testSeq != 'Test733':
                        T2_mask_full[:, :, jj]   = np.flip(np.rot90(aux2_T2_mask_full))
                        T2_mask_single[:, :, jj] = np.flip(np.rot90(aux2_T2_mask_single))
                        T2_mask[:, :,jj]         = T2_mask_single[:, ::-1,jj]
                    else:
                        T2_mask_single[:, :, jj] = (aux2_T2_mask_single)
                        T2_mask[:, :,jj]         = np.rot90(np.rot90(T2_mask_single[::-1, ::-1, jj]))



                    if testSeq == '731' or testSeq == '732':
                        T2_mask_full[:, :, jj]   = aux2_T2_mask_full
                        T2_mask_single[:, :, jj] = np.flip(np.rot90(aux2_T2_mask_single))
                        T2_mask[:, :,jj]         = T2_mask_single[:, ::-1,jj]
                    else:
                        T2_mask_full[:, :, jj]   = aux2_T2_mask_full
                        T2_mask_single[:, :, jj] = aux2_T2_mask_single
                        T2_mask[:, :,jj]         = np.rot90(np.rot90(T2_mask_single[::-1, ::-1, jj]))

                        if subject == 'Subject01':
                            if  testSeq == '733':
                                if dr == 'Day1_Rep1':
                                    T2_mask_single[:, :, jj] = (aux2_T2_mask_single)
                                    T2_mask[:, :, jj]        = np.rot90(np.rot90(T2_mask_single[::-1, ::-1, jj]))
                                else:
                                    T2_mask_single[:, :, jj] = (aux2_T2_mask_single)
                                    T2_mask[:, :, jj]        = np.rot90(T2_mask_single[:, ::-1, jj])
                            else:
                                T2_mask_single[:, :, jj] = (aux2_T2_mask_single)
                                T2_mask[:, :, jj]        = np.rot90(T2_mask_single[:, ::-1, jj])

                        elif subject == 'Subject02':
                            if  testSeq == '733':
                                if dr == 'Day1_Rep1':
                                    T2_mask_single[:, :, jj] = (aux2_T2_mask_single)
                                    T2_mask[:, :, jj]        = np.rot90(T2_mask_single[:, ::-1, jj])
                                else:
                                    T2_mask_single[:, :, jj] = (aux2_T2_mask_single)
                                    T2_mask[:, :, jj]        = np.rot90(T2_mask_single[:, ::-1, jj])


                #plt.figure(figsize=[12, 8])
                ##plt.imshow(mask, cmap='gray')
                #plt.imshow(aux_T2_mask[:, :, 4], cmap='gray')
                #plt.axis('off')

                # Plot per slice
                if plotTest:
                    plt.figure(figsize=[12, 8])
                    for jj in range(nslices):
                        plt.subplot(2,round(nslices/2)+1,jj+1)
                        plt.title(['slice', jj])
                        plt.imshow(T2_mask[:,:,jj], cmap='gray')
                        plt.axis('off')
                        plt.clim(np.min(np.min(imgMSE_final)), np.max(np.max(imgMSE_final)))



                    slice_inf = 0
                    slice_sup = 0
                    sli       = 4
                    fig = plt.figure()
                    # plot_knee_mask = aux_knee_mask.get_fdata()
                    # plot_knee_mask = np.flip(np.rot90(T2_mask))
                    # plot_knee_mask = (np.rot90(T2_mask))
                    # plot_knee_mask = np.flip((T2_mask))
                    # plot_knee_mask = T2_mask
                    for ind in range(slice_inf, slice_sup + 1):
                        # plt.subplot(5, 6, ind + 1)
                        plt.imshow(np.abs(imgMSE_final[:, :, sli,0]), cmap="gist_gray")
                        # plt.imshow(np.rot90(imgMSE_final[:, ::-1, sli,0]), cmap="gist_gray")
                        # plt.imshow(np.rot90(T2_mask[::-1, ::-1, sli]), alpha=.4, cmap="inferno")
                        # plt.imshow(T2_mask[:, ::-1,sli], alpha=.4, cmap="inferno")

                        plt.imshow(T2_mask[:, :,sli], alpha=.4, cmap="inferno")

                    fig = plt.figure()
                    plt.imshow(T2_mask[:, :, 2], alpha=.4, cmap="inferno")

                os.chdir(maskFld)
                # save parameters
                if saveResults:
                    parametersMASK = {'T2_mask_single': T2_mask,'T2_mask_full':T2_mask_full}
                    savemat("roiMSE_py.mat", parametersMASK)

                 # ... 6.2 - segment Knee - Manual Cartilage Segmentation T2 ...
            elif testSegment == 'ManualSegm':

                if maskTest:
                    os.chdir(maskFld)
                    filename = 'Knee_Segment.nii.gz'
                    aux_knee_mask = nib.load(filename)
                    T2_mask = aux_knee_mask.get_fdata()

                    # knee_mask = np.flip(np.rot90(knee_mask[::-1, :, slice_inf]))
                    T2_mask = np.flipud(np.flip(np.rot90(T2_mask[:, :, 0])))

                    hdrMask = aux_knee_mask.header
                    if plotTest:
                        slice_inf = 0
                        slice_sup = 0
                        fig = plt.figure()
                        # plot_knee_mask = aux_knee_mask.get_fdata()
                        # plot_knee_mask = np.flip(np.rot90(plot_knee_mask[:, :, ::-1]))
                        plot_knee_mask = T2_mask
                        for ind in range(slice_inf, slice_sup + 1):
                            # plt.subplot(5, 6, ind + 1)
                            plt.imshow(imgMSE_final[:, :, ind, 1], cmap="gist_gray")
                            # plt.imshow(np.flipud(plot_knee_mask[:, :, ind]), alpha=.4, cmap="inferno")
                            plt.imshow(plot_knee_mask[:, :], alpha=.4, cmap="inferno")

                    os.chdir(maskFld)
                    # save parameters
                    if saveResults:
                        parametersMASK = {'T2_mask': T2_mask}
                        savemat("roiMSE_py.mat", parametersMASK)

                else:
                    os.chdir(maskFld)
                    mask = scipy.io.loadmat('roiMSE_py')
                    T2_mask = mask['T2_mask']

                    if plotTest:
                        fig = plt.figure()
                        plt.imshow(imgMSE_final[:, :, 0, 1], cmap="gist_gray")
                        plt.imshow(T2_mask[:, :], alpha=.4, cmap="inferno")



                # ... 3.1.2 - Number of slices to study ...
            #echo, slice_inf, slice_sup = 1, 1, nslices+1
            #echo, slice_inf, slice_sup = echo-1, slice_inf-1, slice_sup-1


            print('\n\n 6 - Sucessfully finished - inputs of Data (Mask) \n\n')




            # =============================================================================
            #%% --- 6.5 - Read NIFTI data preprocessed
            # =============================================================================

            if niftiTest:


                in_path  = dir_set + '/imagResults_nifti'  # Replace with your input path
                out_path = dir_set + '/imagResults_nifti_Resampled'  # Replace with your input path
                if not os.path.exists(out_path):
                    os.makedirs(out_path)

                # Initialize
                Nx_old = Nx
                Ny_old = Ny
                Nx_new = 512
                Ny_new = 512
                Nx     = Nx_new
                Ny     = Ny_new

                aux_imgMSE_final = np.zeros([Nx,Ny,nslices,nTE])

                for ech in range(0, nTE):
                    data_in_path  = in_path + '/nifti_Data_Knee_ech' + str(ech + 1) + '.nii'
                    data_out_path = out_path + '/nifti_Data_Knee_ech' + str(ech + 1) + '_Resampled.nii'
                    rotated_image = sitk.ReadImage(data_in_path, outputPixelType=sitk.sitkFloat32)

                    # rotated_image = sitk.DICOMOrient(rotated_image, "RAS")
                    resample_img = resample_image(rotated_image, target_size=[str(Nx), str(Ny), '*'])  # Resample to 512x512 inplane resolution, preserving the original number of slices.

                    # Save file in the new folder
                    sitk.WriteImage(resample_img, data_out_path)

                    # Get matrix
                    img_nifti               = nib.load(data_out_path)
                    aux_imgMSE_final[:,:,:,ech] = np.array(img_nifti.dataobj)

                imgMSE_final            = np.rot90(aux_imgMSE_final[:, ::-1, :, :])

                if subject == 'Subject01' or subject == 'Subject02':
                    if testSeq == '733' and dr == 'Day1_Rep1':
                        imgMSE_final = aux_imgMSE_final

                # Figures
                # plt.figure()
                # plt.title("imgMSE from .dat")
                # plt.imshow((imgMSE_final[:, :,0,0]), cmap="gist_gray")
                # plt.show




                # ... 5.3 - Plots
            if plotTest:
                plt.figure(figsize=[12,8])
                plt.title('k-space')
                plt.imshow(abs(aux_k[:,1,:,0,1])**0.2, cmap='gray', origin='lower')
                plt.axis('off')

                # Plot per echo
                plt.figure(figsize=[12, 8])
                for ii in range(nTE):
                    plt.subplot(2,round(nTE/2)+1,ii+1)
                    plt.title(['Echo', ii])
                    plt.imshow(imgMSE_final[:,:,0,ii], cmap='gray')
                    plt.axis('off')
                    plt.clim(np.min(np.min(imgMSE_final)), np.max(np.max(imgMSE_final)))

                # Plot per slice
                plt.figure(figsize=[12, 8])
                for jj in range(nslices):
                    plt.subplot(2,round(nslices/2)+1,jj+1)
                    plt.title(['slice', jj])
                    plt.imshow(imgMSE_final[:,:,jj,0], cmap='gray')
                    plt.axis('off')
                    plt.clim(np.min(np.min(imgMSE_final)), np.max(np.max(imgMSE_final)))

            print('\n\n 4 - Sucessfully finished - load Data MSE & Reconstruction from Nifti \n\n')


            # =============================================================================
            # =============================================================================
            #%% --- 7 - Build dictionary -  T2 estimation
            # =============================================================================
            # =============================================================================

                # ... 7.1 - parameters for dictionary ...
            T1      = 1000                      # ms


            # Get slice profile points
            # Npoints = 50  # Number of points for The SLR profile (lower as possible)
            Npoints = 30    # Number of points for The SLR profile (lower as possible)

            # Performing crop due to noise
            if testCrop:
                crop_ms = 100
                dataVal = 901
                new_nTE = round( crop_ms / ESP)
                if new_nTE < nTE:
                    nTE = new_nTE
                    RFflipAngle = RFflipAngle[0:new_nTE,:]
            else:
                if DenoiseTest:
                    if complexDict:
                        dataVal = 804
                    else:
                        dataVal = 805
                else:
                    if complexDict:
                        dataVal = 803
                    else:
                        dataVal = 801

                # T2 value
            if dataVal == 527 or dataVal == 617 or dataVal == 707 or dataVal == 801 or dataVal == 802:
                aux1_T2 = np.linspace(1,250,250)   # ms
            elif dataVal == 617:  # jun16
                aux1_T2 = np.linspace(1,250,250)   # ms
            elif dataVal == 707:  # jul05
                aux1_T2 = np.linspace(1,250,250)   # ms
            elif dataVal == 801 or dataVal == 901:  # jul09
                aux1_T2 = np.linspace(1,250,250)   # ms
            elif dataVal == 802:
                aux1_T2 = np.linspace(1,250,250)   # ms
            elif dataVal == 803:            # aux1_T2 = np.linspace(1,250,250)   # ms
                aux1_T2 = np.linspace(1,251,251)   # ms
            elif dataVal == 804:            # aux1_T2 = np.linspace(1,250,250)   # ms
                aux1_T2 = np.linspace(1,249,249)   # ms
            elif dataVal == 805:  # aux1_T2 = np.linspace(1,250,250)   # ms
                aux1_T2 = np.linspace(1, 252, 252)  # ms
            T2      = aux1_T2
                        # aux2_T2 = np.linspace(302,800,250)   # ms
            # T2      = np.append(aux1_T2, aux2_T2)

            if dataVal == 527:   # jun27
                B1      = np.linspace(0.65,1.35,71)
            elif dataVal == 617: # jun16
                B1      = np.linspace(0.95,1.25,31)
            elif dataVal == 707: # jul05
                B1      = np.linspace(0.44,0.88,12)
            elif dataVal == 801 or dataVal == 901 or dataVal == 803 or dataVal == 804 or dataVal == 805: # jul09
                B1      = np.linspace(0.90,1.15,26)
            elif dataVal == 802:
                B1      = np.linspace(0.90,1.20,31)

            # Get name to see if dict is needed
            if ConstantFA:
                if SLR_prof:
                    nameMatrix = f"DICT_SliProf_MSE_CRLB_{B1test}_py_ESP_{round(ESP)}ms_nTE_{nTE}_cfA_T2_{T2test}ms_stEx_{round(st_exc * 1e3, 2)}mm_stRF_{round(st_refoc * 1e3, 2)}mm_symm{Tsymm}_sinc{Tsinc}_Npoints{Npoints}_T2min{int(min(T2))}_T2points{T2.shape[0]}_T2max{int(max(T2))}_B1min{int(min(B1) * 100)}_B1points{B1.shape[0]}_B1max{int(max(B1) * 100)}"
                else:
                    nameMatrix = f"DICT_dirac_MSE_CRLB_{B1test}_py_ESP_{round(ESP)}ms_nTE_{nTE}_cfA_T2_{T2test}ms_stEx_{round(st_exc * 1e3, 2)}mm_stRF_{round(st_refoc * 1e3, 2)}mm_symm{Tsymm}_sinc{Tsinc}_Npoints{Npoints}_T2min{int(min(T2))}_T2points{T2.shape[0]}_T2max{int(max(T2))}_B1min{int(min(B1) * 100)}_B1points{B1.shape[0]}_B1max{int(max(B1) * 100)}"
            else:
                if SLR_prof:
                    nameMatrix = f"DICT_SliProf_MSE_CRLB_{B1test}_py_ESP_{round(ESP)}ms_nTE_{nTE}_vfA_T2_{T2test}ms_stEx_{round(st_exc * 1e3, 2)}mm_stRF_{round(st_refoc * 1e3, 2)}mm_symm{Tsymm}_sinc{Tsinc}_Npoints{Npoints}_T2min{int(min(T2))}_T2points{T2.shape[0]}_T2max{int(max(T2))}_B1min{int(min(B1) * 100)}_B1points{B1.shape[0]}_B1max{int(max(B1) * 100)}"
                else:
                    nameMatrix = f"DICT_dirac_MSE_CRLB_{B1test}_py_ESP_{round(ESP)}ms_nTE_{nTE}_vfA_T2_{T2test}ms_stEx_{round(st_exc * 1e3, 2)}mm_stRF_{round(st_refoc * 1e3, 2)}mm_symm{Tsymm}_sinc{Tsinc}_Npoints{Npoints}_T2min{int(min(T2))}_T2points{T2.shape[0]}_T2max{int(max(T2))}_B1min{int(min(B1) * 100)}_B1points{B1.shape[0]}_B1max{int(max(B1) * 100)}"
            nameMatrix = nameMatrix.replace(".", "_")


            dic_exist  = os.path.join(dict_dir, nameMatrix+'.mat')

            if not os.path.exists(dic_exist):
                runDict = 1  # did it run dictionary - '1' yes OR '0' no

                if SLR_prof:
                    FA_refoc = RFflipAngle

                    # ... 7.2 - Get dictionary ...  NOTE: 'FA_exc' in 'rad' & 'FA_refoc' in 'degrees'
                if complexDict:

                    [dict_Phantom, pars] = dict_pars_generator_seq_vJuly25(T1=T1,T2=T2,B1=B1,ESP=ESP,ETL=nTE,
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

                    Dict_Phantom_norm = np.zeros(  (nTE, dict_Phantom.shape[1]) , dtype=np.complex64 )

                else: # with complex numbers
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

                    Dict_Phantom_norm = np.zeros(  (nTE, dict_Phantom.shape[1]) )

                col_T2 = pars[:, 0]  # indexes of T2
                col_B1 = pars[:, 2] * 100  # indexes of B1

                # ... 7.3 - Normalize dictionary ...
                for ii in range(dict_Phantom.shape[1]):
                    Dict_Phantom_norm[:,ii] = dict_Phantom[:,ii]/np.linalg.norm(dict_Phantom[:,ii])


                    # ... 7.4 - Save Dictionary ...
                if saveResults:
                    os.chdir(dict_dir)
                    parametersDICT = {'Dict_Phantom_norm': Dict_Phantom_norm, 'col_T2': col_T2, 'col_B1': col_B1}
                    if ConstantFA:
                        if SLR_prof:
                            nameMatrix = f"DICT_SliProf_MSE_CRLB_{B1test}_py_ESP_{round(ESP)}ms_nTE_{nTE}_cfA_T2_{T2test}ms_stEx_{round(st_exc * 1e3, 2)}mm_stRF_{round(st_refoc * 1e3, 2)}mm_symm{Tsymm}_sinc{Tsinc}_Npoints{Npoints}_T2min{int(min(T2))}_T2points{T2.shape[0]}_T2max{int(max(T2))}_B1min{int(min(B1) * 100)}_B1points{B1.shape[0]}_B1max{int(max(B1) * 100)}"
                        else:
                            nameMatrix = f"DICT_dirac_MSE_CRLB_{B1test}_py_ESP_{round(ESP)}ms_nTE_{nTE}_cfA_T2_{T2test}ms_stEx_{round(st_exc * 1e3, 2)}mm_stRF_{round(st_refoc * 1e3, 2)}mm_symm{Tsymm}_sinc{Tsinc}_Npoints{Npoints}_T2min{int(min(T2))}_T2points{T2.shape[0]}_T2max{int(max(T2))}_B1min{int(min(B1) * 100)}_B1points{B1.shape[0]}_B1max{int(max(B1) * 100)}"
                    else:
                        if SLR_prof:
                            nameMatrix = f"DICT_SliProf_MSE_CRLB_{B1test}_py_ESP_{round(ESP)}ms_nTE_{nTE}_vfA_T2_{T2test}ms_stEx_{round(st_exc * 1e3, 2)}mm_stRF_{round(st_refoc * 1e3, 2)}mm_symm{Tsymm}_sinc{Tsinc}_Npoints{Npoints}_T2min{int(min(T2))}_T2points{T2.shape[0]}_T2max{int(max(T2))}_B1min{int(min(B1) * 100)}_B1points{B1.shape[0]}_B1max{int(max(B1) * 100)}"
                        else:
                            nameMatrix = f"DICT_dirac_MSE_CRLB_{B1test}_py_ESP_{round(ESP)}ms_nTE_{nTE}_vfA_T2_{T2test}ms_stEx_{round(st_exc * 1e3, 2)}mm_stRF_{round(st_refoc * 1e3, 2)}mm_symm{Tsymm}_sinc{Tsinc}_Npoints{Npoints}_T2min{int(min(T2))}_T2points{T2.shape[0]}_T2max{int(max(T2))}_B1min{int(min(B1) * 100)}_B1points{B1.shape[0]}_B1max{int(max(B1) * 100)}"
                    nameMatrix = nameMatrix.replace(".", "_")

                    savemat(nameMatrix+".mat", parametersDICT)

            else:
                runDict = 0   # did it run dictionary - '1' yes OR '0' no

                os.chdir(dict_dir)
                if ConstantFA:
                    if SLR_prof:
                        nameMatrix = f"DICT_SliProf_MSE_CRLB_{B1test}_py_ESP_{round(ESP)}ms_nTE_{nTE}_cfA_T2_{T2test}ms_stEx_{round(st_exc * 1e3, 2)}mm_stRF_{round(st_refoc * 1e3, 2)}mm_symm{Tsymm}_sinc{Tsinc}_Npoints{Npoints}_T2min{int(min(T2))}_T2points{T2.shape[0]}_T2max{int(max(T2))}_B1min{int(min(B1) * 100)}_B1points{B1.shape[0]}_B1max{int(max(B1) * 100)}"
                    else:
                        nameMatrix = f"DICT_dirac_MSE_CRLB_{B1test}_py_ESP_{round(ESP)}ms_nTE_{nTE}_cfA_T2_{T2test}ms_stEx_{round(st_exc * 1e3, 2)}mm_stRF_{round(st_refoc * 1e3, 2)}mm_symm{Tsymm}_sinc{Tsinc}_Npoints{Npoints}_T2min{int(min(T2))}_T2points{T2.shape[0]}_T2max{int(max(T2))}_B1min{int(min(B1) * 100)}_B1points{B1.shape[0]}_B1max{int(max(B1) * 100)}"
                else:
                    if SLR_prof:
                        nameMatrix = f"DICT_SliProf_MSE_CRLB_{B1test}_py_ESP_{round(ESP)}ms_nTE_{nTE}_vfA_T2_{T2test}ms_stEx_{round(st_exc * 1e3, 2)}mm_stRF_{round(st_refoc * 1e3, 2)}mm_symm{Tsymm}_sinc{Tsinc}_Npoints{Npoints}_T2min{int(min(T2))}_T2points{T2.shape[0]}_T2max{int(max(T2))}_B1min{int(min(B1) * 100)}_B1points{B1.shape[0]}_B1max{int(max(B1) * 100)}"
                    else:
                        nameMatrix = f"DICT_dirac_MSE_CRLB_{B1test}_py_ESP_{round(ESP)}ms_nTE_{nTE}_vfA_T2_{T2test}ms_stEx_{round(st_exc * 1e3, 2)}mm_stRF_{round(st_refoc * 1e3, 2)}mm_symm{Tsymm}_sinc{Tsinc}_Npoints{Npoints}_T2min{int(min(T2))}_T2points{T2.shape[0]}_T2max{int(max(T2))}_B1min{int(min(B1) * 100)}_B1points{B1.shape[0]}_B1max{int(max(B1) * 100)}"
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
            # elif accTest == 'noAcc':
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
            slice_inf = 0
            slice_sup = nslices
            #nslices   = slice_sup - slice_inf + 1

            if T2_mask.ndim <= 2: # Number of regions for T2
                nregions = 1
            else:
                nregions = T2_mask.shape[2]

            if complexDict:
                X                = np.zeros((nTE,Nx*Ny,nslices), dtype=np.complex64 )
                T2_dict_map      = np.zeros((Nx,Ny,nslices) , dtype=np.complex64 )
                B1_dict_map      = np.zeros((Nx,Ny,nslices) , dtype=np.complex64 )
                T2_avg_dict      = np.zeros((nslices) , dtype=np.complex64 )
                T2_mask_reshp    = np.zeros((nslices,Nx*Ny) , dtype=np.complex64 )
                avgX             = np.zeros((nTE,nslices) , dtype=np.complex64 )
                StdX             = np.zeros((nTE,nslices) , dtype=np.complex64 )

            else:
                X                = np.zeros((nTE,Nx*Ny,nslices))
                T2_dict_map      = np.zeros((Nx,Ny,nslices))
                B1_dict_map      = np.zeros((Nx,Ny,nslices))
                T2_avg_dict      = np.zeros((nslices))
                T2_mask_reshp    = np.zeros((nslices,Nx*Ny))
                avgX             = np.zeros((nTE,nslices))
                StdX             = np.zeros((nTE,nslices))

                # ... 9.2 - Template_matching of the dictionary ...
            # name_save_T2 = "T2_dictmap_v27_6_py.mat"
            # name_save_B1 = "B1_dictmap_v27_6_py.mat"
            name_save_T2 = "T2_dictmap_t" + str(dataVal) +  "_py.mat"
            name_save_B1 = "B1_dictmap_t" + str(dataVal) +  "_py.mat"

            T2resultsMap_path = os.path.join(dir_Results, name_save_T2)

            if not os.path.exists(T2resultsMap_path):
            # a = 1
            # if a == 1:
                print("Running template matching ...")

                # for aux in range(9, nslices): # cycle in the number of vials
                for nsli in range(0, nslices):  # cycle in the number of vials

                    if nslices == 1:
                        aux_T2_mask = T2_mask
                    else:
                        aux_T2_mask = T2_mask[:, :, nsli]

                    # get reshaped mask
                    T2_mask_reshp[nsli, :] = aux_T2_mask.reshape((1, Nx * Ny), order='C')


                    if testCrop:
                        X_MSE = np.squeeze(imgMSE_final[:, :, nsli, 0:nTE])        # test for only 1 slice
                    else:
                        X_MSE = np.squeeze(imgMSE_final[:, :, nsli, :])            # test for only 1 slice
                    X_MSE = np.transpose(X_MSE, (2, 0, 1))
                    X_MSE = np.reshape(X_MSE, (nTE, Nx*Ny), order='C')

                    if complexDict:
                        X[:,:,nsli] = X_MSE
                    else:
                        X[:,:,nsli] = abs(X_MSE)
                    ind_param   = np.zeros((nslices,X.shape[1]))

                    for i in range(X.shape[1]): #47min
                        print(i)
                        if T2_mask_reshp[nsli,i] == 0:
                            ind_param[nsli,i] = 1
                        else:
                            if complexDict:
                                ind_param[nsli, i] = template_match_vComplex(Dict_Phantom_norm, X[:, i, nsli], col_T2,                                                                             col_B1)
                            else:
                                ind_param[nsli,i]  = template_match_test(  Dict_Phantom_norm,  X[:,i,nsli] , col_T2 , col_B1 )



                    index_Value                     = ind_param[nsli,:]
                    index_Value                     = index_Value.astype(int)

                    if runDict==1: # due to the fact that the vector as different dimentions
                        T2_dict = col_T2[(index_Value)]
                        B1_dict = col_B1[(index_Value)]
                    else:
                        T2_dict = col_T2[(0,index_Value)]
                        B1_dict = col_B1[(0,index_Value)]

                    T2_dict_map [:,:,nsli] = np.reshape(T2_dict, (Nx,Ny), order='C') # [Nx, Ny, nSlices,
                    B1_dict_map [:,:,nsli] = np.reshape(B1_dict, (Nx,Ny), order='C') # [Nx, Ny, nSlices,

                        # ... 9.3. - Save index ...
                if saveResults:
                    os.chdir(dir_Results)
                    parametersINDEX = {'T2_dict_map': T2_dict_map}
                    parametersINDEX_B1 = {'B1_dict_map': B1_dict_map}

                    if niftiTest:  # save results for pre_process
                        savemat("T2_dictmap_py_preProc.mat", parametersINDEX)
                        savemat("B1_dictmap_py_preProc.mat", parametersINDEX_B1)
                    else:
                        savemat(name_save_T2, parametersINDEX)
                        savemat(name_save_B1, parametersINDEX_B1)


            else: # load data
                print("T2_dict_map.mat already exists. Load template matching.")

                os.chdir(dir_Results)
                if niftiTest:  # save results for pre_process
                    Index     = scipy.io.loadmat('T2_dictmap_py_preProc')
                    Index     = scipy.io.loadmat('B1_dictmap_py_preProc')
                else:
                    Index_T2     = scipy.io.loadmat(name_save_T2)
                    Index_B1     = scipy.io.loadmat(name_save_B1)
                T2_dict_map = Index_T2['T2_dict_map']
                B1_dict_map = Index_B1['B1_dict_map']

                # ... 9.3 - Plots ...
            if plotTest:
                fig = plt.figure()
                plt.imshow(T2_dict_map[:, :, 4], cmap="plasma")
                plt.colorbar()
                plt.show
                fig = plt.figure()
                plt.imshow(B1_dict_map[:, :, 4], cmap="plasma")
                plt.colorbar()
                plt.show

            print('\n\n 9 - Sucessfully finished - T2 Dictionary Matching \n\n')



            # =============================================================================
            #%% --- 10 - Figures - Add T2 map on top of the magnitude T2w images
            # =============================================================================

            # 10.1 - Figure for T2
            # ... 10.1 - Figures ...
            if finalPlots:
                fig = plt.figure(62, frameon=False)

                for nsli in range(nslices):

                    plt.subplot(3, round(nslices / 3) + 1, nsli + 1)

                    extent = 0, imgMSE_final.shape[0], 0, imgMSE_final.shape[1]

                    Z1  = abs(imgMSE_final[:, :, nsli, 0]) # 1st echo image
                    im1 = plt.imshow(Z1, cmap=plt.cm.gray, interpolation='nearest',
                                     extent=extent)

                    # Get Value and Standard deviation
                    avgT2vial = np.zeros(nslices)
                    stdT2vial = np.zeros(nslices)

                    # ... 10.2 - Mean & Standart Devixation ...
                    Z2 = abs(T2_dict_map[:, :, nsli])
                    Z2[Z2 == 2] = ['nan']

                    masked_data = np.ma.masked_array(Z2, np.isnan(Z2))
                    avg_Z2 = np.average(masked_data)
                    std_Z2 = np.std(masked_data)
                    T2_Dict_coef = std_Z2 / avg_Z2

                    im2 = plt.imshow(Z2, cmap=plt.cm.jet, interpolation='nearest', extent=extent,vmin=0, vmax=100)
                    try:
                        avgT2vial[nsli] = round(avg_Z2, 2)
                        stdT2vial[nsli] = round(std_Z2, 2)

                    except Exception:
                        # Fallback if rounding fails
                        avgT2vial[nsli] = np.nan  # or set to None, 0, etc.
                        stdT2vial[nsli] = np.nan

                    print("All Cartilage | Slice " + str(nsli) + " | T2avg: " + str(avgT2vial[nsli]) + ", std: " + str(stdT2vial[nsli]))

                    # Figures
                    #
                    #plt.title("Vial T2map - Dictionary, slice: " + str(ind) + ", avg: " + str(round(avg_Z2, 2)) + ", std: " + str(round(std_Z2, 2)))

                    # 10.3 - Save
                    plt.title('Slice ' + str(nsli+1))
                    plt.axis('off')
                    plt.clim(0, 100)
                    plt.show

                    # plt.colorbar()
                    plt.show

                os.chdir(dir_Results)
                if SLR_prof:
                    plt.savefig("Cartilage_SP_T2map_Dictionary_Allslices.png" )
                else:
                    plt.savefig("Cartilage_diracSP_T2map_Dictionary_Allslices.png")




                # 10.2 - Figure for B1
                fig    = plt.figure(64, frameon=False)

                for nsli1 in range(nslices):


                    plt.subplot(3, round(nslices / 3) + 1, nsli1 + 1)

                    extent = 0, imgMSE_final.shape[0], 0, imgMSE_final.shape[1]

                    Z1  = abs(imgMSE_final[:, :, nsli1, 0]) # 1st echo image
                    im1 = plt.imshow(Z1, cmap=plt.cm.gray, interpolation='nearest',
                                     extent=extent)

                    # Get Value and Standard deviation
                    avgB1vial = np.zeros(nslices)
                    stdB1vial = np.zeros(nslices)

                    # ... 10.2 - Mean & Standart Devixation ...
                    Z21 = abs(B1_dict_map[:, :, nsli1])
                    Z21[Z21 == B1[0] * 100] = ['nan']

                    masked_data = np.ma.masked_array(Z21, np.isnan(Z21))
                    avg_Z21 = np.average(masked_data)
                    std_Z21 = np.std(masked_data)
                    B1_Dict_coef = std_Z21 / avg_Z21

                    im21 = plt.imshow(Z21, cmap=plt.cm.jet, interpolation='nearest', extent=extent,vmin=0, vmax=100)
                    try:
                        avgB1vial[nsli1] = round(avg_Z21, 2)
                        stdB1vial[nsli1] = round(std_Z21, 2)

                    except Exception:
                        # Fallback if rounding fails
                        avgB1vial[nsli1] = np.nan  # or set to None, 0, etc.
                        stdB1vial[nsli1] = np.nan

                    print("All Cartilage | Slice " + str(nsli1) + " | B1avg: " + str(avgB1vial[nsli1]) + ", B1 std: " + str(stdB1vial[nsli1]))

                    # Figures
                    #
                    #plt.title("Vial T2map - Dictionary, slice: " + str(ind) + ", avg: " + str(round(avg_Z2, 2)) + ", std: " + str(round(std_Z2, 2)))

                    # 10.3 - Save
                    plt.title('Slice ' + str(nsli1+1))
                    plt.axis('off')
                    plt.clim(90, 140)
                    plt.show

                    plt.colorbar()
                    plt.show

                os.chdir(dir_Results)
                if SLR_prof:
                    plt.savefig("Cartilage_SP_B1map_Dictionary_Allslices.png" )
                else:
                    plt.savefig("Cartilage_diracSP_B1map_Dictionary_Allslices.png")


                print("Data from dir:" + dir_set)

                print('\n\n 10 - Sucessfully finished - Add T2 map on top of the magnitude T2w images \n\n')


                # To create Noise to let me know from the end of processing
                duration = 2000  # milliseconds
                freq = 440  # Hz
                winsound.Beep(freq, duration)

                # pylab.show()

                # print("Template matching complete. Output saved to:", output_file)


