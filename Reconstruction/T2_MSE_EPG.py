#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 28 11:18:25 2021
Tests T2 MSE estimation with EPG for knee protocol
@author: tfernandes

Needs:
 - pip install roipoly
 - pip install sigpy
 - Add toolbox from https://github.com/ut - mri-sim-py in a toolbox folder
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
import nibabel as nib
matplotlib.use('TKAgg')
import matplotlib.pyplot as plt
import ipywebrtc as webrtc
import sys
import math
import time
import statistics
import sigpy
import pydicom
import itk
import tkinter as tk
import winsound

from tkinter import *
from tkinter import filedialog
from roipoly import RoiPoly
from scipy.io import savemat
from statistics import mean
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm
from itkwidgets import view
from IPython.display import display
from tkinter import Tk
from tkinter.filedialog import askdirectory

if PC == 0: # ATTENTION - adapt paths according to your setup
    sys.path.append('Github/Reconstruction/T2_EPG')            # Add functions from mytoolbox
    sys.path.append('Github/Toolboxes/mri-sim-py-master/epg')  # Add toolboxes from https://github.com/utcsilab
    from dict_pars_generator import dict_pars_generator
    from template_match import template_match
    #from epgcpmg import FSE_signal


#%% --- 0.5 - Settings
test        = 1

echo, slice_inf, slice_sup = 1, 1, 30 # Number of slices to study
echo, slice_inf, slice_sup = echo-1, slice_inf-1, slice_sup-1

plotTest      = True     # Get plots 'True' or 'False'
saveResults   = True     # Save results 'True' or 'False'

maskTest      = True     # Get Mask 'True' or 'False'
DictTest      = False    # If already Dictionary 'False' otherwise 'True'
SLR_prof      = True     # Use SLR profile 'True' or 'False'
templateMatch = True     # Do template match 'True'

dataType      = 'dicom'  # Matlab files: 'mat' OR DICOM files: 'dicom'
maskType      = 'NIfTI'  # Nifty: 'NIfTI' OR make mask from data: 'data'


# =============================================================================
# =============================================================================
#%% --- 0 - Select directories - ATTENTION adapth paths according to your setup
# =============================================================================
# =============================================================================

dict_dir = "Github/Data/Dictionaries"
dir_Data = "Github/Examples"

# folder   = "C:/Users/filia/Documents/PhD/Projetos/qMRI/Reconstruction/HSa_REPSTRESS/V002-R Seg OK/Stress/4-Joelho_Dto_Domin"
folder = askdirectory(title='Select Folder') # select directory were folders 'dicom', 'segmentation' & 'results' of data to study are
maskFld  = folder + "/segmentation"

print(folder)

# =============================================================================
# =============================================================================
#%% --- 1 - load Data MSE
# =============================================================================
# =============================================================================

t0 = time.clock()

if dataType == 'mat':
    dir_set = 'Github/Data/template_data' # ATTENTION Change dir according to your setup
    os.chdir(dir_set)

        # ... 1.1 - get file ...
    mat = scipy.io.loadmat('HSA_knee_13032020')

        # ... 1.2 - get data ...
    MSE_TE       = mat['MSE_TE']
    dicom_info   = mat['dicom_info']
    imgMSE_final = mat['imgMSE_final']
    nTE          = int(mat['nechoes'])
    nslices      = int(mat['nslices'])

        # ... 1.3 - set parameters ...
    TE_first    = float(MSE_TE[0][0])  # First TE - in (s)
    ESP         = float(MSE_TE[0][0])  # Spacing Time between Echoes - in (s)
    TE_vector   = np.linspace(TE_first, ESP * nTE, nTE)  # TE vector - in (s)

    nl = imgMSE_final.shape[0]  # Image size

elif dataType == 'dicom':
    # ... 1.1 - get path ...
    # filepath = input("Enter the path of your dicom files: ") # Automatic
    filepath = folder + "/dicom" # dicom file path
    dir_set  = folder + "/results"

        # ... 1.2 - load the data ...
    ind             = 0
    dataset         = []
    aux_data_dicom  = []
    aux_slLoc       = []
    aux_MSE_TE      = []

    for img in os.listdir(filepath):
        # ... 1.2.1 - Get path and dicom files ...
        dicom_file = os.path.join("/", filepath, img)
        dataset.append(pydicom.read_file(dicom_file))

        dicom_info = dataset[ind]
        aux_image_type = dicom_info.ImageType
        image_type = ""
        for ele in aux_image_type:
            image_type += ele

        # ... 1.2.2 - Save Magnitude Images ...
        if image_type == 'ORIGINALPRIMARYM_SEMSE':
            aux_data_dicom.append(dataset[ind].pixel_array)

            aux_slLoc.append(dicom_info.SliceLocation)
            aux_MSE_TE.append(dicom_info.EchoTime)

        print("DICOM reader ite = ", ind+1)
        ind = ind + 1;

    # ... 1.2.3 - Get vectors for slLoc and MSE_TE ...
    slLoc       = np.zeros( (1,len(aux_slLoc)) )
    MSE_TE      = np.zeros( (1,len(aux_MSE_TE)) )
    data_dicom  = np.zeros( (aux_data_dicom[0].shape[0], aux_data_dicom[0].shape[1],len(aux_data_dicom)) )
    ind = 0
    for index in aux_slLoc:
        slLoc[0,ind]  = index
        ind = ind+1
    ind = 0
    for index in aux_MSE_TE:
        MSE_TE[0,ind] = index
        ind = ind+1
    ind = 0
    for index in aux_data_dicom:
        data_dicom[:,:,ind] = index
        ind = ind + 1

    # ... 1.2.4 - Get sequence parameters ...
    slLocs    = np.unique(slLoc)
    nslices   = int(slLocs.shape[0])
    nechoes   = int(dicom_info.EchoTrainLength)
    TR        = int(dicom_info.RepetitionTime)    # Repetition Time in ms
    nTE       = nechoes
    indSl     = np.argsort(slLoc, kind='mergesort')

    imgMSE = data_dicom[:,:,indSl]
    MSE_TE = MSE_TE[0,indSl]

    for ii in range(nslices):
        ind = (nechoes*(ii)+np.linspace(0,nechoes-1,nechoes)).astype(int)
        tmp = imgMSE[:,:,0,ind].astype(int)

        indTE = np.argsort(MSE_TE[0,ind], kind='mergesort')
        imgMSE[:,:,0,ind] = tmp[:,:,indTE].astype(int)

    aux_imgMSE   = np.reshape( imgMSE , (imgMSE.shape[0], imgMSE.shape[1], imgMSE.shape[3]), order="F")
    imgMSE_final = np.reshape( aux_imgMSE , (aux_imgMSE.shape[0], aux_imgMSE.shape[1], nechoes, nslices ), order="F")
    imgMSE_final = np.transpose( imgMSE_final , (0, 1, 3, 2) )

    MSE_TE = np.unique(MSE_TE)

        # ... 1.3 - set parameters ...
    TE_first  = float(MSE_TE[0])  # First TE - in (s)
    ESP       = float(MSE_TE[0])  # Spacing Time between Echoes - in (s)
    TE_vector = np.linspace(TE_first, ESP * nTE, nTE)  # TE vector - in (s)

    nl = imgMSE_final.shape[0]  # Image size

    # ... 1.4 - set Angles ...
FA_exc    = math.pi/2                             # Flip-angle Exciatation - in (rad) - along y
phase_exc = math.pi/2                             # Phase Exciatation - in (rad) - along y

FA_refoc      = 125                                      # Flip-angle Refocosing - in (degrees) - along x
phase_refoc   = np.exp(np.zeros(nTE)/(180/math.pi)*1j)   # Phase Refocosing -  ou zeros(1,nTE);

t1 = time.clock() - t0


fig = plt.figure()
plt.imshow(imgMSE_final[:, :,15,3],
           cmap="gist_gray")  # display image 1st echo random slice with good cartilage display
plt.colorbar()
plt.title("Knee")
plt.show(block=False)

timeVector = np.linspace(0,59,10)
fig = plt.figure()
plt.plot(timeVector,imgMSE_final[187, 134,15,:])
plt.xlabel('Time (ms)', fontsize=25)
plt.ylabel('Norm. Intensity (u.a.)', fontsize=25)
plt.title('Subj5 Right Knee - Signal Decay (154,218), slice 15', fontsize=35)
plt.grid(True)
plt.show



print("Time elapsed to load data: ", t1)
print("1 - Successfully read data")


# =============================================================================
# =============================================================================
#%% --- 2 - inputs of Data (Mask)
# =============================================================================
# =============================================================================

    # ... 2.1 - Create and display images montage ...
# montage_imgMSE = np.array((np.moveaxis(imgMSE_final,-1,2)), dtype=float)     # struct - [Nx, Ny, #echos, #slices]
#
# max_imgMSE = double(max(max(max(abs(montage_imgMSE(:,:, 3,:))))));
# figure;
# montage(montage_imgMSE(:,:, 3,:), [0, max_imgMSE]); % display all slices for 1st echo intensity
#
# montage_imgMSE_disp = permute(imgMSE_final(:,:, slice_inf: slice_sup, 1), [1 2 4 3]); # 1st echo

     # ... 2.2 - segment knee cartilage ...
if maskType == 'data':
    if maskTest: # Only for one Slice TODO - add more slices
        ii = 0
        # knee_mask = np.zeros((imgMSE_final.shape[0],imgMSE_final.shape[1],sli))
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
            knee_mask = roiMSE.get_mask(current_image=imgMSE_final[:,:, ind, 1])
            plt.imshow(knee_mask, alpha=.4, cmap="inferno")
            ii = ii + 1

            # Save image mask
            os.chdir(dir_set)
            plt.savefig("knee_mask_slice_%i.png" % ind)

        os.chdir(dir_set)
        # save parameters
        if saveResults:
            parametersMASK = {'knee_mask': knee_mask}
            savemat("roiMSE_py.mat", parametersMASK)

    else:
        # Plot the slice images if plotTest == 'True'
        if plotTest:
            for ind in range(slice_inf, slice_sup + 1):
                fig = plt.figure(1)
                plt.imshow(imgMSE_final[:, :, ind, echo], cmap="gist_gray")  # display image 1st echo random slice with good cartilage display
                plt.colorbar()
                plt.title("Knee Mask")
                plt.show(block=False)
                os.chdir(dir_set)
                mask = scipy.io.loadmat('roiMSE_py')
                knee_mask = mask['knee_mask']
                plt.imshow(knee_mask, alpha=0.3, cmap="inferno")
                # Save image mask
                os.chdir(dir_set)
                plt.savefig("knee_mask_slice_%i.png" %ind)

        mask = scipy.io.loadmat('roiMSE_py')
        # mask = scipy.io.loadmat('z_knee_mask')
        knee_mask = mask['knee_mask']
        # knee_mask = np.rot90(np.flip(knee_mask, axis=1))
        if plotTest:
            aux_knee_mask = np.rot90(np.flip(knee_mask,axis=1))
            fig = plt.figure(2)
            plt.title("Knee Mask")
            # plt.imshow(aux_knee_mask, cmap="gist_gray")
            plt.imshow(knee_mask, cmap="inferno")
            # plt.imshow(np.reshape(knee_mask, (int(np.sqrt(knee_mask.shape[1])),int(np.sqrt(knee_mask.shape[1])))), cmap="Greys")


if maskType == 'NIfTI':
    if maskTest:
        os.chdir(maskFld)
        filename      = 'Knee_Segment.nii.gz'
        aux_knee_mask = nib.load(filename)
        knee_mask     = aux_knee_mask.get_fdata()
        if slice_inf == slice_sup:
            knee_mask     = np.flip(np.rot90(knee_mask[::-1, :, slice_inf]))
        else:
            knee_mask     = np.flip(np.rot90(knee_mask[::-1, :, slice_inf:slice_sup+1]))

        hdrMask       = aux_knee_mask.header
        if plotTest:
            fig = plt.figure()
            plot_knee_mask = aux_knee_mask.get_fdata()
            plot_knee_mask = np.flip(np.rot90(plot_knee_mask[::-1, :, :]))
            for ind in range(slice_inf, slice_sup + 1):
                plt.subplot(5, 6,ind+1)
                plt.imshow(imgMSE_final[:, :, ind, 1], cmap="gist_gray")
                plt.imshow(plot_knee_mask[:, :, ind], alpha=.4, cmap="inferno")

        os.chdir(dir_set)
        # save parameters
        if saveResults:
            parametersMASK = {'knee_mask': knee_mask}
            savemat("roiMSE_py.mat", parametersMASK)

    else:
        os.chdir(dir_set)
        mask = scipy.io.loadmat('roiMSE_py')
        knee_mask = mask['knee_mask']

        if plotTest:
            fig = plt.figure()
            plt.imshow(imgMSE_final[:, :, 4, 1], cmap="gist_gray")
            plt.imshow(knee_mask[:, :, 4], alpha=.4, cmap="inferno")

print("2 - Successfully read segmented knee")




# save mat files - Data & Mask
#os.chdir(dir_set)
#parametersMAT_imgANDmask = {'imgMSE_final': imgMSE_final, 'knee_mask': knee_mask, 'nslices': nslices, 'nTE': nTE,'ESP':ESP}
#savemat("parametersMAT_imgANDmask.mat", parametersMAT_imgANDmask)


# =============================================================================
# =============================================================================
#%% --- 3 - Build dictionary - bi-exponencial T2 estimation
# =============================================================================
# =============================================================================




# ... 3.1 - parameters for dictionary ...
T1 = 600  # ms
T2 = np.linspace(1, 300, 300)  # ms
B1 = np.linspace(0.6, 1.4, 81)
#T2 = np.linspace(1, 150, 10)  # ms
#B1 = np.linspace(0.6, 1.4, 2)

#If dictTest is True
if DictTest:
    t0 = time.clock()

        # ... 3.2 - Get dictionary .. NOTE: 'FA_exc' in 'rad' & 'FA_refoc' in 'degrees'
    [dict_knee, pars] = dict_pars_generator(T1=T1,T2=T2,B1=B1,ESP=ESP,ETL=nTE,refoc_phase=phase_refoc,phase_exc=np.array([phase_exc]),
                                                   FA_exc=np.array([FA_exc]), FA_refoc=np.array([FA_refoc]), dir_rf=dir_Data, SLR=SLR_prof) # 10min

    col_T2 = pars[:,0]        # indexes of T2
    col_B1 = pars[:,2]*100    # indexes of B1

        # ... 3.3 - Normalize dictionary ...
    Dict_knee_norm = np.zeros((nTE,dict_knee.shape[1]))

    for ii in range(dict_knee.shape[1]):
        Dict_knee_norm[:,ii] = dict_knee[:,ii]/np.linalg.norm(dict_knee[:,ii])

    if plotTest:
        fig = plt.figure()
        plt.plot(abs(Dict_knee_norm[:,2]))

    t1 = time.clock() - t0
    print("Time elapsed per slice: ", t1)

        # ... 3.4 - Save Dictionary ...
    if saveResults:
        os.chdir(dict_dir)
        parametersDICT = {'Dict_knee_norm': Dict_knee_norm,'col_T2': col_T2, 'col_B1': col_B1}
        nameMatrix = f"DICT_MSE_py_ESP_{round(ESP)}ms_nTE_{nTE}_fA_{FA_refoc}"
        savemat(nameMatrix+".mat", parametersDICT)


else:
    os.chdir(dict_dir)
    nameMatrix = f"DICT_MSE_py_ESP_{round(ESP)}ms_nTE_{nTE}_fA_{FA_refoc}"
    Dict = scipy.io.loadmat(nameMatrix)
    Dict_knee_norm  = Dict['Dict_knee_norm'] # TODO alterar retirando o shortTE
    col_T2          = Dict['col_T2']
    col_B1          = Dict['col_B1']


print("3 - Successfully got Dictionary")


timeVector = np.linspace(0,59,10)
fig = plt.figure()
plt.plot(timeVector,Dict_knee_norm[:, :])
plt.xlabel('Time (ms)', fontsize=25)
plt.ylabel('Norm. Intensity (u.a.)', fontsize=25)
plt.title('Dictionary', fontsize=35)
plt.grid(True)
plt.show





# =============================================================================
# =============================================================================
#%% --- 4 - T2 Dictionary Matching
# =============================================================================
# =============================================================================



# ... 4.1 - initialize parameters ...
nslices         = slice_sup - slice_inf + 1
X               = np.zeros((nTE,nl*nl,nslices))
T2_dict_map     = np.zeros((nl,nl,nslices))
ind_param       = np.zeros((nslices, int(nl*nl)))

itt_slice = 0

# plt.figure()
# plt.plot(knee_mask_reshp)
# plt.show()

    # ... 4.2 - Template_matching of the dictionary ...
if templateMatch:

    # =============================================================================
    # =============================================================================

    # aux_kMask = knee_mask[:, :, 2]
    # knee_mask_reshp = aux_kMask.reshape((1, nl * nl), order='C')
    #
    # X_MSE     = np.squeeze(imgMSE_final[282-1, 263-1, 2, :])
    # X_test    = abs(X_MSE)
    # ind_param = template_match(Dict_knee_norm, X_test)
    #
    # fig = plt.figure()
    # plt.plot(timeVector,X_test/ np.linalg.norm(X_test), label='Data')
    # plt.plot(timeVector,Dict_knee_norm[:, ind_param], label='Best Dict Match')
    # plt.xlabel('Time (ms)', fontsize=25)
    # plt.ylabel('Norm. Intensity (u.a.)', fontsize=25)
    # plt.title('Dictionary fit', fontsize=35)
    # plt.grid(True)
    # plt.legend(['Data', 'Best Dict Match'])
    # plt.show

    # =============================================================================
    # =============================================================================

    for ind in range(slice_inf,slice_sup+1,1):
        if slice_inf == slice_sup:
            knee_mask_reshp = knee_mask.reshape((1, nl * nl), order='C')
        else:
            aux_kMask       = knee_mask[:, :, itt_slice]
            knee_mask_reshp = aux_kMask.reshape((1, nl * nl), order='C')

        X_MSE = np.squeeze(imgMSE_final[:,:,ind,:])
        X_MSE = np.transpose(X_MSE, (2, 0, 1))
        X_MSE = np.reshape(X_MSE, (nTE, nl*nl), order='C')
        X[:,:,itt_slice] = abs(X_MSE)

        for i in range(X.shape[1]): #47min - run over points within the mask for all echoes, and checks best fit with dictionary
            print(i)
            if knee_mask_reshp[0,i] == 0:
                ind_param[itt_slice,i] = 1
            else:
                ind_param[itt_slice,i] = template_match(  Dict_knee_norm,  X[:,i,itt_slice]  )

        index_Value                 = ind_param[itt_slice,:]
        index_Value                 = index_Value.astype(int)
        T2_dict                     = col_T2[(0,index_Value-1)]+1
        T2_dict_map [:,:,itt_slice] = np.reshape(T2_dict, (nl,nl), order='C') # [Nx, Ny, nSlices,
        itt_slice = itt_slice + 1

    # ... 4.3. - Save index ...
    if saveResults:
        os.chdir(dir_set)
        parametersINDEX = {'T2_dict_map': T2_dict_map}
        savemat("T2_dictmap_py.mat", parametersINDEX)

else:
    os.chdir(dir_set)
    Index     = scipy.io.loadmat('T2_dictmap_py')
    T2_dict_map = Index['T2_dict_map']

    # ... 4.3 - Plots ...
if plotTest:
    fig = plt.figure()
    plt.imshow(T2_dict_map[:, :, 4], cmap="plasma")
    plt.colorbar()
    plt.show

print("4 - Successfully performed Template Matching")




# =============================================================================
# =============================================================================
#%% --- 5 - Figures - Add T2 map on top of the magnitude T2w images
# =============================================================================
# =============================================================================

itt_slice=0
#Z2_3D = np.zeros((nl,nl,nslices))
Z2_3D = np.zeros((nl,nl,slice_sup-slice_inf+1))
Z2_3D_nifti = np.zeros((nl,nl,slice_sup-slice_inf+1))

for ind in range(slice_inf, slice_sup + 1):
#for ind in range(27):
    fig = plt.figure(ind+61, frameon=False)

    extent = 0, imgMSE_final.shape[0], 0, imgMSE_final.shape[1]

    Z1  = imgMSE_final[:, :, ind, echo]
    im1 = plt.imshow(Z1, cmap=plt.cm.gray, interpolation='nearest',
                     extent=extent)

    Z2 = T2_dict_map[:, :, itt_slice]
    Z2[Z2==2] = ['nan']

        # ... 5.2 - Mean & Standart Devixation ...
    if sum(sum(np.isin(Z2,T2))) > 0:    # checks if there is mask
        masked_data = np.ma.masked_array(Z2, np.isnan(Z2))
        avg_Z2 = np.average(masked_data)
        std_Z2 = np.std(masked_data)
        T2_Dict_coef = std_Z2 / avg_Z2

        im2 = plt.imshow(Z2, cmap=plt.cm.jet, interpolation='nearest',
                         extent=extent)
        plt.title("Cartilage T2map - Dictionary, slice: " +str(ind+1) + ", avg: "+str(round(avg_Z2, 2)) + ", std: " + str(round(std_Z2, 2)) )
        plt.colorbar()
        plt.clim(0, nTE * ESP)
        plt.show
    else:
        plt.title("Cartilage T2map - Dictionary, slice: " +str(ind+1) + ", no mask for knee cartilage")
        plt.show

    # Save for 3D representation
    Z2_3D[:,:,ind] = Z2

    if saveResults:
        os.chdir(dir_set)
        slicPlot = ind +1
        plt.savefig("Cartilage_T2map_Dictionary_slice_%i.png" %slicPlot)

        if maskType == 'NIfTI':
            T2mapImage = nib.Nifti1Image(np.flip(np.rot90(Z2[::-1, :])), affine=np.eye(4))
            T2mapText  = ("Cartilage_T2map_Dictionary_slice_%i.nii.gz" %slicPlot)
            nib.save(T2mapImage, os.path.join(dir_set, T2mapText))
            Z2_3D_nifti[:, :, ind] = np.flip(np.rot90(Z2[::-1, :]))

    itt_slice = itt_slice + 1

if saveResults and maskType == 'NIfTI':
    os.chdir(dir_set)
    T2mapALL_image = nib.Nifti1Image(Z2_3D_nifti, affine=np.eye(4))
    T2mapALL_text  = ("Cartilage_ALLT2map_Dictionary.nii.gz")
    nib.save(T2mapALL_image, os.path.join(dir_set, T2mapALL_text))

    # no inversion of image
#    T2mapALL_image_test = nib.Nifti1Image(Z2_3D, affine=np.eye(4))
#    T2mapALL_text_test  = ("Cartilage_ALLT2map_Dictionary_test.nii.gz")
#    nib.save(T2mapALL_image_test, os.path.join(dir_set, T2mapALL_text_test))

    # =============================================================================
    # =============================================================================
    # %% --- 6 -  3D MAP: Interactive rendering of 𝑇2 maps
    # =============================================================================
    # =============================================================================

if saveResults: # read data with 'MSE_recondata_3D_T2maps.m'
    parMAT = {'data':Z2_3D}
    savemat("test_mat_Maps.mat", parMAT)

print(folder)
print("------- CODE RUN SUCCESSFULLY --------------")


duration = 1000  # milliseconds
freq = 440  # Hz
winsound.Beep(freq, duration)