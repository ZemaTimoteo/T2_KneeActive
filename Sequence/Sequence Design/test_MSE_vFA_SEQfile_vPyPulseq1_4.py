#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May 5 12:33:25 2020
@author: tfernandes @ISR-IST - Uni Lisbon
"""

# =============================================================================

# %% --- 0 - Import functions
import os
import sys
import scipy.io
import numpy as np
import matplotlib.pyplot as plt
import math
import warnings
import scipy.io
import time
from scipy.io import savemat
from scipy.interpolate import interp1d
from pathlib import Path

import pypulseq

from pypulseq.Sequence.sequence import Sequence
from pypulseq.calc_rf_center import calc_rf_center
from pypulseq.make_adc import make_adc
from pypulseq.make_delay import make_delay
from pypulseq.make_extended_trapezoid import make_extended_trapezoid
from pypulseq.make_sinc_pulse import make_sinc_pulse
from pypulseq.make_gauss_pulse import make_gauss_pulse
from pypulseq.make_trapezoid import make_trapezoid
from pypulseq.opts import Opts
from pypulseq.calc_duration import calc_duration
from pypulseq.make_block_pulse import make_block_pulse


# %% ===========================================================================
user_profile = os.environ.get("USERPROFILE") or os.environ.get("HOME")  # Check user home directory
if user_profile:
    # For Windows, Documents is under the user's profile folder.
    documents_folder = Path(user_profile) / "Documents"

#####
# %% --- 0 - Settings ---
#####

save_flag   = True  # save file

selectDir   = ''     # selected directory for save

# Types of Test
diffST_Test      = True      # True - Different slice thickness for exc and refoc | False - Same slice thickness for exc and refoc
SincRFTest       = True      # True - Sinc | False - Gauss
symmetricRF      = False     # Test for symmetric RF pulse 'True' OR 'False' - asymmetric
intercalateLines = True      # True - Intercalate lines to avoid Magnetization Transfer | False - not intercalated lines
unfoldTest       = False     # Increase the Nx in readout to avoid the Folding
sliceGAP         = True      # True - interval between Slices | False - Slices are contiguous

# Sequence Acceleration
accTest = 'LORAKS'   # 'SENSE' OR 'GRAPPA' OR 'LORAKS' OR 'noAcc'

# Other Parameters
plotTest    = False  # Plot sequence
reportTest  = False  # Report the sequence
rf_pulses   = False  # Only perform rf_pulses
ktrajTest   = False  # Test K_trajectory


# Type of Test
kneeTest    = True  # False - Phantom | True - Knee

# --------------------------------------------------

if intercalateLines:
    broTest = True  # 'True' - Bit-reversed slack OR 'False' - Alternated order
else:
    broTest = False  # 'True' - Bit-reversed slack OR 'False' - Alternated order

#####
# --- 1 - Set sequence limits ---
#####

# optimization Parameters - input parameters
experiment_id = 1000
n_slices      = 15      # slices
T2            = 45     # in ms


# Optimized vFA
if T2==8: # for Target T2=8ms
    # with Slice Profile
    FA_best = np.array([162.515,      140.8719,      138.9243,      135.6338,      135.0291,      130.3957,
                        135.61,      149.3292,      151.4102,      23.68045,       36.5044,      55.59507,
                        78.71435,       97.3574,      108.3908,      109.5766,      97.33363,      83.30018,
                        77.32394,       86.5669])

    n_echo  = FA_best.shape[0]             # #echoes
    TE      = 8.5e-3                       # in s (min: 8e-3) (min: 8.95e-3 for Nx = 256)
    TR      = 2608e-3                      # in s


elif T2==45: # for Target T2=45ms
    FA_best = np.array([150.1127,      120.8003,      114.9025,      112.5717,      112.9336,       111.505,      110.0765,
                        108.6479,      107.2193,      109.3717,      111.5241,      113.6766,       115.829,      116.1909,
                        112.9718,      111.5433,      116.3814,      115.8481,      121.5815,      147.0103,      168.8582])
    # TE      = 8.5e-3                        # in s (min: 8e-3) (min: 8.95e-3 for Nx = 256)
    TE      = 8.5e-3                        # in s (min: 8e-3) (min: 8.95e-3 for Nx = 256)
    TR      = 3073e-3                       # in s - Original da optimizaçao
    n_echo  = FA_best.shape[0]              # #echoes

# rf_flip_Vector = FA_best*math.pi/180
rf_flip_Vector = FA_best


# RF angles
TBWP_ex      = 4             # Time Bandwidth Product in RF exc
TBWP_ref     = 2             # Time Bandwidth Product in RF refoc
rf_ex_angle  = 90            # RF excitation - in degrees

# scan parameters
maxGrad = 32    # in mT/m (max - 32)
maxSlew = 140   # in T/m/s (max - 140)

dG      = 160e-6           # rise time (s)


# Slices
#n_slices        = 2     # #slices
slice_thickness = 3e-3  # in m - (min: 2.6e-3)
slicevectReal   = np.linspace(1,n_slices,n_slices)

if diffST_Test:
    slice_thickness_exc   = round(slice_thickness,4)
    slice_thickness_refoc = round(slice_thickness * 3,4)
else:
    slice_thickness_exc   = slice_thickness
    slice_thickness_refoc = slice_thickness

pe_type = 'linear'

i_raster_time = 100000
DT = 1 / i_raster_time  # in s

if kneeTest: # Geometry for Knee - image parameters
    fov     = 170e-3  # in m
    Nx      = 230
else: # Geometry for Phantom
    fov     = 256e-3  # in m
    # Nx      = 400
    Nx      = 256

Ny      = Nx

rf_ex_freq_offset  = 0   # Frequency | offset in Hertz(Hz)
rf_ref_freq_offset = 0   # Frequency | offset in Hertz(Hz)

delta_k = 1 / fov
k_width = Nx * delta_k

gamma = 42.54e6  # Hz/Tesla

if accTest == 'SENSE':
    R   = 2         # Acceleration Factor for
    Ny_real = Ny/R  # real Ny for sequence

elif accTest == 'GRAPPA': # +25 & -25 to centre k-space, and the rest (200-50=150) equaly sampled
    R               = 2                         # Acceleration Factor
    fullLin         = 24
    percFullKspace  = int(fullLin/Nx)           # Percentage of Full K-space

elif accTest == 'LORAKS':
    R               = 8                         # Acceleration Factor - normalmente 4
    percFullKspace  = 0.17                      # Percentage of Full K-space
    fullLin         = int(Nx*percFullKspace)    # Fully Sampled K-space center with x% of Nx

elif accTest == 'noAcc':
    R               = 1                         # Acceleration Factor - normalmente 4
    percFullKspace  = 1                         # Percentage of Full K-space
    fullLin         = int(Nx*percFullKspace)    # Fully Sampled K-space center with x% of Nx


if symmetricRF:
    symmetricVal = 0.5
else:
    symmetricVal = 0.85  # the bigger then 0.5 the more to the right, the less the 0.5 the more to the left

if sliceGAP:
    sliceGAPsiz = 3 * 1e-3  # slice Gap, in (m)
else:
    sliceGAPsiz = 0          # slice Gap, in (m)

# timing
adc_dead_time    = 10e-6                                                    # ADC dead Time
rf_dead_time     = 100e-6                                                   # RF dead Time
rf_ringdown_time = 20e-6                                                    # RF Ringdown Time
delay_adc        = 20e-6                                                    # Delay ADC

# Bandwidth and Dwell Time
rBW              = 100000                                                    # Bandwidth in Hz
# rBW              = 75000                                                    # Bandwidth in Hz
dwell_time       = np.round(1/rBW,decimals=7)
# dwell_time       = 20e-6                                                    # Dwell Time (Time between ADC Points)

# readout_time     = 2.5e-3 + 2 * adc_dead_time                             # time of readout
if unfoldTest == True:  # double the Aquisition FOV in readout Direction
    origNx      = Nx
    Nx          = Nx * 2
    dwell_time  = dwell_time/2                                              # Dwell Time (Time between ADC Points)
elif unfoldTest == False:
    origNx      = Nx

readout_time     =  ( math.ceil(dwell_time * i_raster_time * Nx) + math.ceil(2*delay_adc * i_raster_time)) /  i_raster_time   # time of readout (in s)

# readout_time     = Nx*system
t_ex             = 2.8e-3                                                   # time of excitation rf_pulse   - 2.9376e-3 (in s)
t_exwd           = t_ex + rf_ringdown_time + rf_dead_time                   # time of excitation window| time of excitation + ringdown time + dead time
t_ref            = 1.6e-3                                                   # time of refocusing rf_pulse  - 1.5232e-3
t_refwd          = t_ref + rf_ringdown_time + rf_dead_time                  # time of refocusing window| time of refocusing + ringdown time + dead time
t_sp             = 0.5 * (TE - readout_time - t_refwd)                      # time of spoiler| TE>dG & cannot be smaller than readout_time + refwindow
t_spex           = 0.5 * (TE - t_refwd) - (1-symmetricVal) * t_ex - 0.5 * (
                        rf_ringdown_time + rf_dead_time)                    # time of spoiler excitation | it is with this complex calcule because of what is time of the rfexcitation
                                                                            # (moved t_exwd out of the parentisis - Original expression: 0.5 * (TE - t_exwd - t_refwd),
                                                                            # have to be careful with deadtime and ringtime )

if TE < 4 * dG + readout_time + t_refwd:
    raise ValueError('TE is too small compared with rise time (dG)')


############################################################################################################################################
################################################### no input after this point ##############################################################
############################################################################################################################################


# timing in number of points
n_adc_dead_time    = math.ceil(adc_dead_time * i_raster_time)
n_rf_dead_time     = math.ceil(rf_dead_time * i_raster_time)
n_rf_ringdown_time = math.ceil(rf_ringdown_time * i_raster_time)
n_delay_adc        = math.ceil(delay_adc * i_raster_time)

n_readout_time     = math.ceil(readout_time * i_raster_time)
n_t_ex             = math.ceil(t_ex * i_raster_time)
n_t_exwd           = math.ceil(t_exwd * i_raster_time)
n_t_ref            = math.ceil(t_ref * i_raster_time)
n_t_refwd          = math.ceil(t_refwd * i_raster_time)

n_t_sp             = math.ceil(t_sp * i_raster_time)
n_t_spex           = math.ceil(t_spex * i_raster_time)

n_dG               = math.ceil(dG * i_raster_time)
n_TE               = math.ceil(TE * i_raster_time)


#####
# --- 2 - Instantiation and gradient limits ---
#####
system = Opts(max_grad=maxGrad, grad_unit='mT/m', max_slew=maxSlew, slew_unit='T/m/s', rf_ringdown_time=n_rf_ringdown_time / i_raster_time,
              rf_dead_time=n_rf_dead_time / i_raster_time, adc_dead_time=n_adc_dead_time / i_raster_time)
seq    = Sequence(system)

# I need to check the manually inputed inverse raster time and the actual value are the same.
assert 1 / i_raster_time == seq.grad_raster_time, "Manualy inputed inverse raster time does not match the actual value."



fsp_r = 1
fsp_s = 0.5



#####
# %% --- 3 - RF pulses
#####
# Excitation - CPMG
rf_ex_phase  = np.pi / 2
rf_ref_phase = 0

        # 3.1 - RF excitation
flip_ex = rf_ex_angle * np.pi / 180

rf_ex, gz, gz_reph = make_sinc_pulse(flip_angle=flip_ex, system=system, duration=n_t_ex / i_raster_time,
                               slice_thickness=slice_thickness_exc,
                               center_pos=symmetricVal, apodization=0.5, time_bw_product=TBWP_ex,
                               freq_offset=rf_ex_freq_offset, phase_offset=rf_ex_phase,
                               return_gz=True)


if calc_duration(rf_ex) - t_exwd > 1e-5:  # allow for time check of rf_ex & slice selection flat time
    t_exwd = calc_duration(rf_ex)
    n_t_exwd = math.ceil(t_exwd * i_raster_time)

if save_flag:
    subfolder1 = 'T2KneeActive'
    subfolder2 = 'rf_pulses'
    dir_rf = Path(documents_folder) / subfolder1 / subfolder2

    # Create the subfolders if they don't exist
    try:
        dir_rf.mkdir(parents=True, exist_ok=True)  # 'parents=True' creates intermediate directories if needed
    except Exception as e:
        print(f"Error creating subfolders: {e}")

    # save
    os.chdir(dir_rf)
    parametersMATexc = {'rf_ex': rf_ex, 't_ex': t_ex, 'G_ex': gz.amplitude}
    nameRFexcsave    = "rf_pulses_exc{}".format(rf_ex_angle) + "_symm{}".format(symmetricVal) + "_st{}".format(slice_thickness_exc*1e3) + "mm"
    nameRFexcsave    = nameRFexcsave.replace(".", "_")
    nameRFexcsave    = nameRFexcsave + ".mat"
    savemat(nameRFexcsave, parametersMATexc)

        # 3.2 - RF refocusing
flip_ref = np.zeros([1, n_echo], dtype=np.complex64)
rf_ref   = np.empty(n_echo, dtype=object)
gz_ref   = np.empty(n_echo, dtype=object)

for k_echo in range(n_echo):  # cicle per TR each with one phase (line of k_space) and multiple Echos
    flip_ref[0,k_echo] = rf_flip_Vector[k_echo] * np.pi / 180
    if SincRFTest:
        testSinc = 1
        aux_rf_ref, aux_gz_ref , _  = make_sinc_pulse(flip_angle=abs(flip_ref[0,k_echo]), system=system, duration=n_t_ref / i_raster_time,
                                        slice_thickness=slice_thickness_refoc,
                                        apodization=0.5, time_bw_product=TBWP_ref,
                                        freq_offset=rf_ref_freq_offset,phase_offset=rf_ref_phase, use='refocusing',
                                        return_gz=True)
    else:
        testSinc = 0
        aux_rf_ref, aux_gz_ref, _   = make_gauss_pulse(flip_angle=abs(flip_ref[0,k_echo]), system=system, duration=n_t_ref / i_raster_time,
                                             slice_thickness=slice_thickness_refoc,
                                             apodization=0.5, time_bw_product=TBWP_ref,
                                             freq_offset=rf_ref_freq_offset,phase_offset=rf_ref_phase,
                                             use='refocusing',
                                             return_gz=True)

    if calc_duration(aux_rf_ref) - t_refwd > 1e-5:  # allow for time check of rf_ref & slice selection flat time - Allways the same time
        t_refwd   = calc_duration(aux_rf_ref)
        n_t_refwd = math.ceil(t_refwd * i_raster_time)

    rf_ref[k_echo] = aux_rf_ref
    gz_ref[k_echo] = aux_gz_ref

    # Save rf pulses - save '.mat'
    if save_flag:
        os.chdir(dir_rf)
        parametersMATref = {'rf_ref': rf_ref[k_echo], 't_ref': t_ref, 'G_ref': gz_ref[k_echo].amplitude}  # original até sept24
        # parametersMATref = {'rf_ref': rf_ref[k_echo], 't_ref': t_ref, 'G_ref': gz.amplitude}                # New due to line 311
        # nameRFrefsave    = "rf_pulses_ref{}".format(round(rf_flip_Vector[k_echo])) + "_{}.mat".format(int(math.modf(rf_flip_Vector[k_echo])[0]*100))

        nameRFrefsave = "rf_pulses_ref" + str(rf_flip_Vector[k_echo])
        if nameRFrefsave[-2:] == '_0':
            nameRFrefsave = nameRFrefsave[:-2]
        nameRFrefsave = nameRFrefsave + "_sinc{}".format(testSinc) + "_st{}".format(slice_thickness_refoc*1e3) + "mm"
        nameRFrefsave = nameRFrefsave.replace(".", "_")
        nameRFrefsave = nameRFrefsave + ".mat"
        savemat(nameRFrefsave, parametersMATref)

print("\n\n ----- rf_pulses saved -----  ")

if rf_pulses:
    print("\n\n quit")
    quit()



#####
# %% --- 4 - Gradients Slice Selection
#####
gs_ex  = make_trapezoid(channel='z', system=system, amplitude=gz.amplitude, flat_time=n_t_exwd / i_raster_time,
                        rise_time=n_dG / i_raster_time)
gs_ref = make_trapezoid(channel='z', system=system, amplitude=gz_ref[0].amplitude, flat_time=n_t_refwd / i_raster_time,
                        rise_time=n_dG / i_raster_time)

# Slice Selection spoiling / considering possible if the case of asymmetric RFexc
if symmetricRF:
    ags_ex  = gs_ex.area / 2
else:
    gs_ex_rise_area = (gs_ex.area-gs_ex.flat_area)/2
    ags_ex  = gs_ex.flat_area*(1-symmetricVal) + gs_ex_rise_area

# ags_ex  = - gs_ex.area / 2 # has to be negative

amplitudeTest = (ags_ex * (1 + fsp_s)) / (n_t_spex / i_raster_time - n_dG / i_raster_time)
# maxGrad*42.54*1e3> amplitudeTest

gs_spr  = make_trapezoid(channel='z', system=system, area=ags_ex * (1 + fsp_s), duration=n_t_sp / i_raster_time,
                         rise_time=n_dG / i_raster_time)
gs_spex = make_trapezoid(channel='z', system=system, area=ags_ex * fsp_s,       duration=n_t_spex / i_raster_time,
                         rise_time=n_dG / i_raster_time)


print("\n\n ... gz rise_time =", np.round(gz.rise_time*1e3,3), " ms | Amplitude =", np.round(gz.amplitude/gamma,3)*1e3, "mT/m | slew rate = ", np.round(gz.amplitude/gamma/gz.rise_time,3), " mT/m/ms")

#####
# %% --- 5 - ADCs / Readouts
#####
# Readout gradient and ADC
gr_acq             = make_trapezoid(channel='x', system=system, flat_area=k_width, flat_time=n_readout_time / i_raster_time, rise_time=n_dG / i_raster_time)  # Gradient Readout
n_gr_acq_flat_time = math.ceil(gr_acq.flat_time * i_raster_time)

adc                = make_adc(num_samples=Nx, duration=(n_gr_acq_flat_time - 2 * n_delay_adc)/ i_raster_time, delay=n_delay_adc / i_raster_time)

adc.dwell = np.round(adc.dwell,decimals=7)
# adc.dwell = math.ceil(adc.dwell*i_raster_time) / i_raster_time
if adc.dwell != dwell_time:
    raise ValueError('Dwell Time is not respect due to Nx or readout_time')
    
#####
# %% --- 6 - Spoilers - Gradient Readout
#####
# RO spoiling
gr_spr  = make_trapezoid(channel='x', system=system, area=gr_acq.area * fsp_r, duration=n_t_sp / i_raster_time, rise_time=n_dG / i_raster_time)
gr_spex = make_trapezoid(channel='x', system=system, area=gr_acq.area * (1 + fsp_r), duration=n_t_spex / i_raster_time, rise_time=n_dG / i_raster_time)

agr_spr = gr_spr.area

# Prephasing gradient in RO direction
agr_preph = gr_acq.area / 2 + agr_spr

gr_preph  = make_trapezoid(channel='x', system=system, area=agr_preph, duration=n_t_spex / i_raster_time, rise_time=n_dG / i_raster_time)


#####
# %% --- 7 - Defining Phase Areas
#####
# n_ex = math.floor(Ny / n_echo)
# pe_steps = np.arange(1, n_echo * n_ex + 1) - 0.5 * n_echo * n_ex - 1
# if divmod(n_echo, 2)[1] == 0:
#     pe_steps = np.roll(pe_steps, -round(n_ex / 2))
# pe_order = pe_steps.reshape((n_ex, n_echo), order='F').T
# phase_areas = pe_order * delta_k

# Phase encoding lobe calculations
pe_steps = np.floor(np.arange(1, Ny + 1) - 0.5 * Ny - 1)
phase_areas = pe_steps * delta_k
n_phA = len(phase_areas)


#####
# %% --- 8 - Split gradients and recombine into blocks - Gradient_Readout & Gradient_Slice_Selection
#####

# ----------------  Gz ----------------------------------------------
# gs1 : ramp up of gs_ex
gs1_times   = [0, gs_ex.rise_time]
gs1_amp     = [0, gs_ex.amplitude]
n_gs1_times = [0, math.ceil(gs_ex.rise_time * i_raster_time)]
n_gs1_amp   = [0, math.ceil(gs_ex.amplitude * i_raster_time)]
gs1         = make_extended_trapezoid(channel='z', times=[n_gs1_times[0] / i_raster_time, n_gs1_times[1] / i_raster_time],
                                      amplitudes=[n_gs1_amp[0] / i_raster_time, n_gs1_amp[1] / i_raster_time])

area_gs1    = gs_ex.amplitude*gs_ex.rise_time/2

# gs2 : flat part of gs_ex
gs2_times = [0, gs_ex.flat_time]
gs2_amp   = [gs_ex.amplitude, gs_ex.amplitude]
gs2       = make_extended_trapezoid(channel='z', times=gs2_times, amplitudes=gs2_amp)

area_gs2  = gs_ex.amplitude*gs_ex.flat_time

# gs3 : Bridged slice pre-spoiler
gs3_times = [0, gs_spex.rise_time, gs_spex.rise_time + gs_spex.flat_time,
             gs_spex.rise_time + gs_spex.flat_time + gs_spex.fall_time]
gs3_amp   = [gs_ex.amplitude, gs_spex.amplitude, gs_spex.amplitude, gs_ref.amplitude]
gs3       = make_extended_trapezoid(channel='z', times=gs3_times, amplitudes=gs3_amp)

area_gs3  =  ( gs_spex.rise_time*(gs_ex.amplitude - gs_spex.amplitude)/2 + gs_spex.rise_time*gs_spex.amplitude ) + gs_spex.amplitude * gs_spex.flat_time + ( (gs_ref.amplitude-gs_spex.amplitude)*gs_spex.fall_time/2 + gs_spex.amplitude * gs_spex.fall_time )

# gs4 : Flat slice selector for pi-pulse
gs4_times = [0, gs_ref.flat_time]
gs4_amp   = [gs_ref.amplitude, gs_ref.amplitude]
gs4       = make_extended_trapezoid(channel='z', times=gs4_times, amplitudes=gs4_amp)
calc_duration(gs4)

area_gs4  = gs_ref.flat_time*gs_ref.amplitude

# gs5 : Bridged slice post-spoiler
gs5_times = [0, gs_spr.rise_time, gs_spr.rise_time + gs_spr.flat_time,
             gs_spr.rise_time + gs_spr.flat_time + gs_spr.fall_time]
gs5_amp   = [gs_ref.amplitude, gs_spr.amplitude, gs_spr.amplitude, 0]
gs5       = make_extended_trapezoid(channel='z', times=gs5_times, amplitudes=gs5_amp)
calc_duration(gs5)

area_gs5  =  ( gs_spr.rise_time*(gs_ref.amplitude - gs_spr.amplitude)/2 + gs_spr.rise_time*gs_spr.amplitude) + gs_spr.flat_time * gs_spr.amplitude + gs_spr.amplitude*gs_spr.fall_time/2


# gs7 : The gs3 for next pi-pulse
gs7_times = [0, gs_spr.rise_time, gs_spr.rise_time + gs_spr.flat_time,
             gs_spr.rise_time + gs_spr.flat_time + gs_spr.fall_time]
gs7_amp   = [0, gs_spr.amplitude, gs_spr.amplitude, gs_ref.amplitude]
gs7       = make_extended_trapezoid(channel='z', times=gs7_times, amplitudes=gs7_amp)

area_gs7  =  gs_spr.rise_time*gs_spr.amplitude/2 + gs_spr.flat_time * gs_spr.amplitude + ( (gs_ref.amplitude-gs_spr.amplitude)*gs_spr.fall_time/2 + gs_spr.amplitude*gs_spr.fall_time )

# calculate moment for rf exc and rf refoc
test_ZeroMoment_90_180_ADC = (area_gs2 * (1 - symmetricVal) ) + area_gs3 + area_gs4/2 - area_gs4/2 - area_gs5
test_ZeroMoment_90_180_ADC == 0

# calculate moment for rf refoc and rf refoc
test_ZeroMoment_ADC_180_ADC = area_gs7 + area_gs4/2 - area_gs4/2 - area_gs5
test_ZeroMoment_ADC_180_ADC == 0


# ----------------  Gz ----------------------------------------------
# ----------------  Gx ----------------------------------------------


# gr3 : pre-readout gradient
gr3      = gr_preph
area_gr3 = gr_preph.area

# gr5 : Readout post-spoiler
gr5_times = [0, gr_spr.rise_time, gr_spr.rise_time + gr_spr.flat_time,
             gr_spr.rise_time + gr_spr.flat_time + gr_spr.fall_time]
gr5_amp   = [0, gr_spr.amplitude, gr_spr.amplitude, gr_acq.amplitude]
gr5       = make_extended_trapezoid(channel='x', times=gr5_times, amplitudes=gr5_amp)
area_gr5  = gr_spr.rise_time*gr_spr.amplitude/2 + gr_spr.flat_time * gr_spr.amplitude + ( gr_spr.fall_time*(gr_spr.amplitude-gr_acq.amplitude)/2 + gr_spr.fall_time*gr_acq.amplitude)

# gr6 : Flat readout gradient
gr6_times = [0, readout_time]
gr6_amp   = [gr_acq.amplitude, gr_acq.amplitude]
gr6       = make_extended_trapezoid(channel='x', times=gr6_times, amplitudes=gr6_amp)
area_gr6  = readout_time*gr_acq.amplitude

# gr7 : the gr3 for next repeat
gr7_times = [0, gr_spr.rise_time, gr_spr.rise_time + gr_spr.flat_time,
             gr_spr.rise_time + gr_spr.flat_time + gr_spr.fall_time]
gr7_amp   = [gr_acq.amplitude, gr_spr.amplitude, gr_spr.amplitude, 0]
gr7       = make_extended_trapezoid(channel='x', times=gr7_times, amplitudes=gr7_amp)
area_gr7  = ( gr_spr.rise_time*(gr_spr.amplitude-gr_acq.amplitude)/2 + gr_spr.fall_time*gr_acq.amplitude) + gr_spr.flat_time * gr_spr.amplitude + gr_spr.fall_time*gr_spr.amplitude/2

# calculate moment for rf exc and rf refoc
test_ZeroMoment_Gx_90_180_ADC = area_gr3 - area_gr5 - area_gr6/2
test_ZeroMoment_Gx_90_180_ADC == 0


# ----------------  Gx ----------------------------------------------



#####
#%% --- 9 - Slice Ordering
#####

if intercalateLines:
    if broTest: # Bit-reversed ordering (BRO) slack
        n_sli_vect = []
        bro_vector = np.array([1, 9, 5, 13, 3, 11, 7, 15, 2, 10, 6, 14, 4, 12, 8, 16])  # Slack of slices | Bit-reversed Ordering (bro)
        if n_slices/(max(bro_vector)+1) > 1:
            for i in range(0,1+math.floor(n_slices/(max(bro_vector)+1))):
                n_sli_vect = np.append(n_sli_vect,(bro_vector+bro_vector.shape[0]*(i) ))
                bro_vector = np.flip(bro_vector)
        else:
            n_sli_vect = bro_vector

    else:      # Alternated order
        a_list = list(range(1, n_slices+1))
        aux_a = a_list[::2]
        if len(a_list)>1:
            aux_b = a_list[1::2]
        else:
            aux_b = []
        n_sli_vect = aux_a + aux_b

else:
    n_sli_vect = np.linspace(1, n_slices, n_slices)


#####
#%% --- 10 - Calculate timing and delays
#####
# delay_TR : delay at the end of each MSE pulse train (i.e. each TR)
total_time_ex  = gs1.tt[-1] + gs2.tt[-1] + gs3.tt[-1]
total_time_ref = gs4.tt[-1] + gs5.tt[-1] + readout_time + gs7.tt[-1]
total_time_end = gs4.tt[-1] + gs5.tt[-1]
TE_train       = total_time_ex + n_echo*total_time_ref + total_time_end    # in s

if broTest: # add slack to TR
    n_slicesTR    = n_sli_vect.shape[0]
    delay_noSLice = make_delay(TE_train)
else:
    n_slicesTR = n_slices

TR_fill  = (TR - n_slicesTR * TE_train)   # time in (s)
TR_fill  = system.grad_raster_time * round(TR_fill / system.grad_raster_time)

if TR_fill < 0:
    TR_fill = 1e-3
    warnings.warn(f'TR too short, adapted to include all slices to: {1e3 * n_slicesTR * (TE_train + TR_fill)} ms')
else:
    print(f'time p/ ETL used = {int(1000 * TE_train)} ms of TR = {int(1e3 * TR)} ms, Time available p/TR -> TR - (ETL time + waiting for other slices) -> = {int(1e3 * TR_fill)} ms')

# Delay
delay_TR = make_delay(TR_fill)

# set times in pointsF
n_TE_train = math.ceil(TE_train * i_raster_time)
n_TR_fill  = math.ceil(TR_fill * i_raster_time)
n_delay_TR = math.ceil(delay_TR.delay * i_raster_time)



#####
#%% --- 11 - Parallel Imaging w/ SENSE factor R
#####

if accTest == 'SENSE':
    n_phA = np.round(int(n_phA / R))
    phase_areas = phase_areas[::R]

elif accTest == 'GRAPPA':
    phase_areas = pe_steps * delta_k
    n_phA = len(phase_areas)
    restLin = (origNx - fullLin) / R
    aux_centerKspac = fullLin / 2
    iniFullcent = int(n_phA / 2 - aux_centerKspac - 1)  # Ini of Fully Sampled
    endFullcent = int(n_phA / 2 + aux_centerKspac)  # End of Fully Sampled

    n_phA = int(fullLin + ((origNx - fullLin) / R))
    Ny_real = n_phA  # real Ny for sequence
    aux_pA_1 = phase_areas[0:iniFullcent:R]  # Get phases of AF for 1st interval
    iniKfull = int(aux_pA_1.shape[0])  # Point where it start fully sampled k-space
    aux_pA_2 = phase_areas[iniFullcent + 1:endFullcent + 1]  # Get phases of AF for 2nd interval
    endKfull = iniKfull + int(aux_pA_2.shape[0])  # Point where it ends fully sampled k-space
    aux_pA_3 = phase_areas[endFullcent + 2:int(phase_areas.shape[0]):R]  # Get phases of AF for 3rd interval

    phase_areas = np.concatenate((aux_pA_1, aux_pA_2, aux_pA_3), axis=0)

elif accTest == 'LORAKS':
    # AuxVariables
    n_phA_original = len(phase_areas)
    aux_centerKspac = fullLin / 2
    iniKfull = int((n_phA / 2 - 1) - aux_centerKspac)  # Ini of Fully Sampled
    endKfull = int((n_phA / 2 - 1) + aux_centerKspac)  # End of Fully Sampled

    Half_Ny = int(math.floor((Ny - fullLin) / 2) - 1)
    newNy = int(Half_Ny / R)
    mask_Kspace = np.zeros((Ny, n_echo))
    phase_areas_matrix = np.zeros((newNy * 2 + fullLin, n_echo))

    # Center fully sampled
    kspaceSampl_2 = np.linspace(iniKfull, endKfull - 1,
                                (endKfull - 1) - (iniKfull) + 1)  # Get phases of AF for 2nd interval (Full K-space)

    for echo in range(n_echo):  # cicle per TR each with one phase (line of k_space) and multiple Echos
        aux_kspaceSampl_1 = np.random.permutation(Half_Ny)  # Permute in the first third of the kspace
        aux_kspaceSampl_3_a = np.linspace(endKfull, Ny - 1, (Ny - 1) - (endKfull) + 1)  # Get final third in points
        aux_kspaceSampl_3_b = np.random.permutation(Half_Ny)  # Permute in the final third of kspace

        kspaceSampl_1 = np.sort(aux_kspaceSampl_1[0:newNy])
        kspaceSampl_3 = aux_kspaceSampl_3_a[np.sort(aux_kspaceSampl_3_b[0:newNy])]
        kspaceSample_all = np.concatenate(
            (kspaceSampl_1.astype(int), kspaceSampl_2.astype(int), kspaceSampl_3.astype(int)), axis=0)

        mask_Kspace[kspaceSample_all, echo] = 1  # Mask for recon with LORAKS
        phase_areas_matrix[:, echo] = phase_areas[kspaceSample_all]  # Matrix with Phases

    n_phA = len(phase_areas_matrix)
    Ny_real = n_phA


#####
# %% --- 12 - Define sequence blocks/Readouts + Create '.seq' file
#####

# SAR_check = True
#
# while SAR_check:
for ph_A in range(n_phA):  # cicle per TR each with one phase (line of k_space) and multiple Echos
    for s in range(0, len(n_sli_vect)): # For each slice (interleaved)
        if n_sli_vect[s] in slicevectReal:

            # rf_ex.freq_offset  = gs_ex.amplitude * ( slice_thickness + sliceGAPsiz ) * (n_sli_vect[s] - (n_slices - 1) / 2)                         # frequency offset for rf pulse - excitation
            rf_ex.freq_offset  = gs_ex.amplitude * ( slice_thickness + sliceGAPsiz ) * ((n_sli_vect[s]-1) - (n_slices - 1) / 2)                       # frequency offset for rf pulse - excitation
            rf_ex.phase_offset = rf_ex_phase - 2 * np.pi * rf_ex.freq_offset * calc_rf_center(rf_ex)[0]                                               # Phase offset for rf pulse - excitation

            seq.add_block(gs1)
            seq.add_block(gs2, rf_ex)
            #seq.add_block(gs2)
            seq.add_block(gs3, gr3)

            for k_echo in range(n_echo):  # For each TR

                #rf_ref[k_echo].freq_offset  = gs_ref.amplitude * ( slice_thickness + sliceGAPsiz ) * (n_sli_vect[s] - (n_slices - 1) / 2)    # frequency offset for rf pulse - reference
                rf_ref[k_echo].freq_offset  = gs_ref.amplitude * ( slice_thickness + sliceGAPsiz ) * ( (n_sli_vect[s]-1) - (n_slices - 1) / 2)    # frequency offset for rf pulse - reference
                rf_ref[k_echo].phase_offset = rf_ref_phase - 2 * np.pi * rf_ref[k_echo].freq_offset * calc_rf_center(rf_ref[k_echo])[0]          # Phase offset for rf pulse - reference

                if accTest == 'LORAKS':
                    phase_area = phase_areas_matrix[ph_A, k_echo]
                else:
                    phase_area = phase_areas[ph_A]

                # Make phase encoding gradients (they are reused for all readouts)
                gp_pre = make_trapezoid(channel='y', system=system, area=phase_area,  duration=n_t_sp / i_raster_time, rise_time=n_dG / i_raster_time)
                gp_rew = make_trapezoid(channel='y', system=system, area=-phase_area, duration=n_t_sp / i_raster_time, rise_time=n_dG / i_raster_time)
                seq.add_block(gs4, rf_ref[k_echo])
                seq.add_block(gs5, gr5, gp_pre)

                seq.add_block(gr6, adc)

                seq.add_block(gs7, gr7, gp_rew)

            # Spoillers
            seq.add_block(gs4)
            seq.add_block(gs5)

        else: # Get delay for when no slice is going to be measured
            seq.add_block(delay_noSLice)
        print( f'-> Slice = {int(s+1)} / {int(len(n_sli_vect))} | in Ny_pos = {int(ph_A+1)} / {int(n_phA)} ...')

    seq.add_block(delay_TR)




# Duration of sequence in ms
seqDur = seq.duration()
print("\n\n .. seq duration .. ", seqDur[0], "in s")
print(f'time p/ ETL used = {int(1000 * TE_train)} ms of TR = {int(1e3 * TR)} ms, Time available p/TR -> TR - (ETL time + waiting for other slices) -> = {int(1e3 * TR_fill)} ms')

#####
# %% --- 13 - B1+rms - control for not exceeding the allowed radiation for the System
#####

# ----------
# Pulseq Implementation from Matlab
dt     = 1e-06   # (in s)
rf_rms = 0

# resample the pulse to a resonable time array
nn_ex  = round(t_ex/dt)
t      = ( np.linspace(0,(nn_ex-1),nn_ex) + 0.5 ) * dt    # in (s)

rfs    = np.interp(t,rf_ex.t,rf_ex.signal)               # in (Hz)
rfs_sq = np.multiply(rfs,np.conj(rfs))

total_energy  = sum(rfs_sq)*dt                             # uT s
peak_pwr      = max(rfs_sq)
aux_rf_rms_ex = np.sqrt(total_energy/t_ex)              # uT^1/2

rf_rms = ( ( aux_rf_rms_ex**2 * t_ex ) )  # uT s

nn_ref     = round(t_ref / dt)
t_intr_ref = (np.linspace(0, (nn_ref - 1), nn_ref) + 0.5) * dt

for k_echo in range(n_echo):  # cicle per TR each with one phase (line of k_space) and multiple Echos
    rfs_ref    = np.interp(t_intr_ref, rf_ref[k_echo].t, rf_ref[k_echo].signal)             # (in Hz)
    rfs_ref_sq = np.multiply(rfs_ref, np.conj(rfs_ref))

    total_energy_ref = sum(rfs_ref_sq) * dt   # uT s
    total_energy     = total_energy + total_energy_ref
    peak_pwr_ref     = max(rfs_ref_sq)
    aux_rf_rms_ref   = np.sqrt(total_energy_ref / t_ref)

    aux_rf_rms =  ( ( aux_rf_rms_ref ** 2 * t_ref) )   # uT s
    rf_rms     = rf_rms + aux_rf_rms                   # uT s

# B1_rms+ for one ETL
aux_rf_rms_perTR = np.sqrt(rf_rms / TE_train)   # Hz
aux_rf_rms = np.sqrt( (rf_rms *  n_slices * n_phA )/ seqDur[0])   # Hz

# Get for all slices and real Ny
B1_rms_perTR = aux_rf_rms_perTR  / gamma  *1e6  # uT s (gamma - Hz/Tesla)
B1_rms = aux_rf_rms  / gamma  *1e6  # uT s (gamma - Hz/Tesla)


print("\n\n .. B1 + rms (from Pulseq)... ", np.round(B1_rms,3), "in uT s")


# ----------
# My  Implementation
B1plus_rf_ex  = np.max(rf_ex.signal) / gamma * 1e6  # units (uT)
b1plus_t_ex   = flip_ex / (
                 np.max(rf_ex.signal) * 2 * np.pi) * 1e3                # time for specific area:  gamma*2pi*B1*time = flip - units (ms) | gamma from the expression cuts with the gamma needed to divided from the signal
sumb1plus_rms = (B1plus_rf_ex ** 2) * b1plus_t_ex * n_slices * n_phA    # (in uT^2 ms) use of n_phA and not Ny because is the true RFexc and RFrefoc needed
sumb1plus_rms_test = (B1plus_rf_ex ** 2) * t_ex * n_slices * n_phA      # (in uT^2 s) use of n_phA and not Ny because is the true RFexc and RFrefoc needed

for k_echo in range(n_echo):  # cicle per TR each with one phase (line of k_space) and multiple Echos
    B1plus_rf_ref  = np.max(rf_ref[k_echo].signal) / gamma * 1e6                                         # units (uT)
    # print("\n\n .. max rfpower of echo ",k_echo, " value:", np.round(B1plus_rf_ref, 3), "in uT")
    b1plus_t_refoc = flip_ref[0,k_echo] / (np.max(rf_ref[k_echo].signal) * 2 * np.pi) * 1e3              # units (ms)
    sumb1plus_rms  = sumb1plus_rms   +  (  (B1plus_rf_ref ** 2) * b1plus_t_refoc  ) * n_slices * n_phA   # (in uT^2 ms) use of n_phA and not Ny because is the true RFexc and RFrefoc needed

    sumb1plus_rms_test  = sumb1plus_rms_test   +  (  (B1plus_rf_ref ** 2) * t_ref  ) * n_slices * n_phA  # (in uT^2 ms) use of n_phA and not Ny because is the true RFexc and RFrefoc needed


b1Plus_rms = np.abs(np.sqrt(sumb1plus_rms / (seqDur[0] * 1e3)))      # units (uT)

#b1Plus_rms_test = np.abs(np.sqrt(sumb1plus_rms_test / (seqDur[0])))  # units (uT)
#b1Plus_rms      = np.abs(np.sqrt(sumb1plus_rms / (10 * 1e3)))     # units (uT)
#b1Plus_rms_test = np.abs(np.sqrt(sumb1plus_rms_test / 10))        #

# print("\n\n .. B1 + rms .. ", np.round(b1Plus_rms,3), "in uT")

# if b1Plus_rms > 5:
#     raise ValueError('B1+rms is higher than threshold - Increase TR')



#####
# %% --- 14 - Plots & Report & K-traject
#####
if plotTest:
    seq.plot()
    # seq.plot(time_range=[0, 0.02 * TR])


if reportTest:
    print(seq.test_report())
    print(seq.check_timing())    

if ktrajTest:
    # Calculate trajectory
    t0 = time.clock()

    ktraj_adc, ktraj, t_excitation, t_refocusing, t_adc = seq.calculate_kspace()
    t1 = time.clock() - t0
    print("Time elapsed per ktraj calculate: ", t1)
    # Plot k-spaces
    time_axis = np.arange(1, ktraj.shape[1] + 1) * system.grad_raster_time
    # plt.figure()
    # plt.plot(time_axis, ktraj.T)  # Plot entire k-space trajectory
    # plt.plot(t_adc, ktraj_adc[0], '.')  # Plot sampling points on kx-axis
    plt.figure()
    plt.plot(ktraj[0], ktraj[1], 'b', ktraj_adc[0], ktraj_adc[1], 'r.')  # 2D plot
    plt.axis('equal')
    plt.show()

#####
# %% --- 15 - save file '.seq'
#####

#print("\n\n .. gz rise_time =", np.round(gz.rise_time*1e3,3), " ms | Amplitude =", np.round(gz.amplitude/gamma,3)*1e3, "mT/m | slew rate = ", np.round(gz.amplitude/gamma/gz.rise_time,3), " mT/m/ms")
#print("\n\n .. gz Ref rise_time =", np.round(aux_gz_ref.rise_time*1e3,3), " ms | Amplitude ref =", np.round(aux_gz_ref.amplitude/gamma,3)*1e3, "mT/m | slew rate ref = ", np.round(aux_gz_ref.amplitude/gamma/aux_gz_ref.rise_time,3), " mT/m/ms")


if save_flag:
    dir = selectDir
    os.chdir(selectDir)

    folder_name = 'test%s_MSE_vFA_T2_%sms' % (experiment_id, T2) + '_acc_' + accTest
    seq_name = 'test%s_pypul_1_4_2_MSE_vFA_T2-%sms_nEchos-%s_nSlices-%s_TE-%sms_TR-%sms_FOV-%smm_Nx-%s_Ny-%s_gm-%s_sm-%s_sTexc-%smm_sTrefoc-%smm' % (
    experiment_id, T2, n_echo, n_slices, round(TE * 1e3), round(TR * 1e3), round(fov * 1e3), Nx, Ny, maxGrad,maxSlew,
    round(slice_thickness_exc * 1e3), round(slice_thickness_refoc * 1e3))
    seq_name = seq_name + '_acc_' + accTest

    if os.path.isdir(folder_name):
        print("\n\n - Directory already exist - ")
    else:
        os.makedirs(folder_name)
    os.chdir(folder_name)

    # get definitions
    # seq.set_definition('FOV', [fov, fov, slice_thickness])
    # seq.set_definition('Name', 'MSE_test')
    seq.write(seq_name)

    # Add definitions
    # seq.definitions()

    # save '.mat' file with parameters
    if accTest == 'SENSE':        # for Parallel Imaging - SENSE
        parametersMAT = {'test': experiment_id, 'TR': TR, 'TE': TE, 'DT': DT, 'nslices': n_slices,
                         'st_exc': slice_thickness_exc, 'st_refoc': slice_thickness_refoc,
                         'max_grad': maxGrad, 'max_slew': maxSlew, 'flipAngle': rf_ex_angle, 'rf_flipAngle': rf_flip_Vector,
                         'nEcohs': n_echo, 'Nx': Nx, 'Ny': Ny, 'Ny_real': Ny_real, 'FOV': fov, 'duration': seqDur[0],
                         'testSinc':testSinc,'symmetricVal':symmetricVal,'sliceGAP':sliceGAP,
                         'accTest': accTest,'n_sli_vect':n_sli_vect,'slicevectReal':slicevectReal,
                         'testIntercalate':intercalateLines,'b1Plus_rms':B1_rms,'broTest': broTest,
                         'delta_k': delta_k, 'k_width': k_width, 'R': R,'T2test':T2,
                         'sliceGAPsiz':sliceGAPsiz, 'unfoldTest':unfoldTest, 'dwellTime': dwell_time, 'rBW': rBW}
    elif accTest == 'GRAPPA':  # for Parallel Imaging - GRAPPA
        parametersMAT = {'test': experiment_id, 'TR': TR, 'TE': TE, 'DT': DT, 'nslices': n_slices,
                         'st_exc': slice_thickness_exc, 'st_refoc': slice_thickness_refoc,
                         'max_grad': maxGrad, 'max_slew': maxSlew, 'flipAngle': rf_ex_angle, 'rf_flipAngle': rf_flip_Vector,
                         'nEcohs': n_echo, 'Nx': Nx, 'Ny': Ny, 'Ny_real': Ny_real, 'FOV': fov,
                         'duration': seqDur[0], 'delta_k': delta_k, 'k_width': k_width, 'R': R,
                         'testSinc':testSinc,'symmetricVal':symmetricVal,'sliceGAP':sliceGAP,
                         'accTest': accTest,'n_sli_vect':n_sli_vect,'slicevectReal':slicevectReal,
                         'testIntercalate':intercalateLines,'b1Plus_rms':B1_rms,'broTest': broTest,
                         'fullLin': fullLin, 'iniKfull': iniKfull, 'endKfull': endKfull,'T2test':T2,
                         'sliceGAPsiz':sliceGAPsiz, 'unfoldTest':unfoldTest, 'dwellTime': dwell_time, 'rBW': rBW}
    elif accTest == 'LORAKS':  # for Sequence Acceleration - LORAKS
        parametersMAT = {'test': experiment_id, 'TR': TR, 'TE': TE, 'DT': DT, 'nslices': n_slices,
                         'st_exc': slice_thickness_exc, 'st_refoc': slice_thickness_refoc,
                         'max_grad': maxGrad, 'max_slew': maxSlew, 'flipAngle': rf_ex_angle, 'rf_flipAngle': rf_flip_Vector,
                         'nEcohs': n_echo, 'Nx': Nx, 'Ny': Ny, 'Ny_real': Ny_real, 'FOV': fov,
                         'duration': seqDur[0], 'delta_k': delta_k, 'k_width': k_width, 'R': R,
                         'testSinc':testSinc,'symmetricVal':symmetricVal,'sliceGAP':sliceGAP,
                         'accTest': accTest,'n_sli_vect':n_sli_vect,'slicevectReal':slicevectReal,
                         'testIntercalate':intercalateLines,'b1Plus_rms':B1_rms,
                         'mask_Kspace':mask_Kspace,'percFullKspace':percFullKspace,'broTest': broTest,
                         'fullLin': fullLin, 'iniKfull':iniKfull+1,'endKfull':endKfull,'T2test':T2,
                         'sliceGAPsiz':sliceGAPsiz, 'unfoldTest':unfoldTest, 'dwellTime': dwell_time, 'rBW': rBW}
    else:  # for normal acquisition
        parametersMAT = {'test': experiment_id, 'TR': TR, 'TE': TE, 'DT': DT, 'nslices': n_slices,
                         'st_exc': slice_thickness_exc, 'st_refoc': slice_thickness_refoc,
                         'max_grad': maxGrad, 'max_slew': maxSlew, 'flipAngle': rf_ex_angle, 'rf_flipAngle': rf_flip_Vector,
                         'nEcohs': n_echo, 'Nx': Nx, 'Ny': Ny, 'FOV': fov, 'duration': seqDur[0], 'delta_k': delta_k,
                         'testSinc':testSinc,'symmetricVal':symmetricVal,'sliceGAP':sliceGAP,
                         'accTest': accTest,'n_sli_vect':n_sli_vect,'slicevectReal':slicevectReal,
                         'testIntercalate':intercalateLines,'b1Plus_rms':B1_rms,'broTest': broTest,
                         'k_width': k_width, 'T2test':T2,'sliceGAPsiz':sliceGAPsiz, 'unfoldTest':unfoldTest, 'dwellTime': dwell_time, 'rBW': rBW}
    savemat("sequence_info.mat", parametersMAT)


    print("\n\n ----- Seq file saved -----  ", seq_name)


# seq.test_report()
a=2