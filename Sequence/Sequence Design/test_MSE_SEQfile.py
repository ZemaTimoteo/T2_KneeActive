#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May 5 12:33:25 2020
@author: tfernandes @ISR-IST - Uni Lisbon
"""

# =============================================================================

#%% --- 0 - Import functions
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

PC = 0      # Seia=1 or myPC=0

if PC == 0:
    import pypulseq

    from pypulseq.Sequence.sequence import Sequence
    from pypulseq.calc_rf_center import calc_rf_center
    from pypulseq.make_adc import make_adc
    from pypulseq.make_delay import make_delay
    from pypulseq.make_extended_trapezoid import make_extended_trapezoid
    from pypulseq.make_sinc_pulse import make_sinc_pulse
    from pypulseq.make_gauss_pulse import make_gauss_pulse
    from pypulseq.make_trap_pulse import make_trapezoid
    from pypulseq.opts import Opts
    from pypulseq.calc_duration import calc_duration
    from pypulseq.make_block_pulse import make_block_pulse

    sys.path.append('Github/Toolboxes/py2jemris')

#from seq2xml import seq2xml

#%% ===========================================================================



#####
#%% --- 0 - Settings ---
#####
test       = 33

save_flag  = False    # save file

plotTest   = True    # Plot sequence
reportTest = False   # Report the sequence
ConvertXML = False   # Convert seq2xml
rf_pulses  = False   # Only perform rf_pulses
ktrajTest  = False   # Test K_trajectory

# Signal Pseudo-steady-state | Flip Back Pulse
FBP_test   = False

# Parallel Imaging
SENSE      = False   # Apply SENSE
GRAPPA     = True    # Apply GRAPPA


#####
# --- 1 - Set sequence limits ---
#####

# optimization Parameters
TE       = 8e-3        # in s (min: 8e-3) (min: 8.95e-3 for Nx = 256)
n_echo   = 10          # #echoes
rf_flip  = 151         # RF refocusing - in degrees
TR       = 2367e-3     # in s


# scan parameters
maxGrad  = 32                # in mT/m (max - 32) | 32
maxSlew  = 130               # in T/m/s (max - 140) | 120
#dG       = 100e-6           # rise time (s)
dG_preph = 200e-6            # rise time for preph (s)
dG       = 180e-6            # rise time (s) - only this one 10 / 2 /2022 - change in gr_preph dG


n_slices = 1                # #slices

rf_ex_angle  = 90            # RF excitation - in degrees
rf_flip_N_1  = 170           # RF refocusing - in degrees - of N + 1 (for FPB)
rf_ref_angle = rf_flip

T1kneeMax = 1100e-3 #(2056e-3)/5          # Max T1 for Cartilage
# TR        = T1kneeMax*5    # in s
TE_eff    = 60e-3            # Effective TE in s
pe_type   = 'linear'

i_raster_time = 100000
DT = 1/i_raster_time         # in s

# Geometry - image parameters
fov             = 256e-3     # in m
Nx              = 256
Ny              = 256
slice_thickness = 5e-3     # in m - (min: 2.6e-3)

# aux_calc
k0      = round(TE_eff / TE)
delta_k = 1 / fov
k_width = Nx * delta_k

gamma   = 42.54e6 # Hz/Tesla

if SENSE:
    R   = 2         # Acceleration Factor for
    Ny_real = Ny/R  # real Ny for sequence

if GRAPPA: # +25 & -25 to centre k-space, and the rest (200-50=150) equaly sampled
    R       = 2         # Acceleration Factor
    fullLin = 24


#####
# --- 2 - Instantiation and gradient limits ---
#####

system = Opts(max_grad=maxGrad, grad_unit='mT/m', max_slew=maxSlew, slew_unit='T/m/s', rf_ringdown_time=100e-6,
              rf_dead_time=100e-6, adc_dead_time=10e-6)
seq = Sequence(system)

# I need to check the manually inputed inverse raster time and the actual value are the same.
assert 1 / i_raster_time == seq.grad_raster_time, "Manualy inputed inverse raster time does not match the actual value."


if isinstance(rf_flip, int):
    rf_flip = np.zeros(n_echo) + rf_flip

# timing
readout_time = 2.5e-3 + 2 * system.adc_dead_time                        # time of readout
# t_ex         = 2.5e-3                                                 # time of excitation rf_pulse   - 2.9376e-3
t_ex        = 1.7e-3                                                    # time of excitation rf_pulse   - 2.9376e-3
t_exwd       = t_ex + system.rf_ringdown_time + system.rf_dead_time     # time of excitation window| time of excitation + ringdown time + dead time
# t_ref        = 1.7e-3                                                 # time of refocusing rf_pulse  - 1.5232e-3
t_ref        = 1.25e-3                                                  # time of refocusing rf_pulse  - 1.5232e-3
t_refwd      = t_ref + system.rf_ringdown_time + system.rf_dead_time    # time of refocusing window| time of refocusing + ringdown time + dead time
t_spex       = 0.5 * (TE - t_exwd - t_refwd)                            # time of spoiler excitation| - gs3
t_sp         = 0.5 * (TE - readout_time - t_refwd)                      # time of spoiler| TE>dG & cannot be smaller than readout_time + refwindow - gs5

if TE - (0.5*t_exwd + t_refwd + readout_time*0.5) - (t_sp+t_spex) > 0.000000001:
    raise ValueError('TE does not match the time for pulses and spoilers')

if TE < 4*dG + readout_time + t_refwd:
    raise ValueError('TE is too small compared with rise time (dG)')


# timing in number of points
n_TE           = math.ceil(TE * i_raster_time)
n_readout_time = math.ceil(readout_time * i_raster_time)
n_t_ex         = math.ceil(t_ex * i_raster_time)
n_t_exwd       = math.ceil(t_exwd * i_raster_time)
n_t_ref        = math.ceil(t_ref * i_raster_time)
n_t_refwd      = math.ceil(t_refwd * i_raster_time)
n_t_sp         = math.ceil(t_sp * i_raster_time)
n_t_spex       = math.ceil(t_spex * i_raster_time)
n_dG           = math.ceil(dG * i_raster_time)
n_dG_preph     = math.ceil(dG_preph * i_raster_time)


fsp_r = 1
fsp_s = 0.5



#####
#%% --- 3 - RF pulses
#####
rf_ex_phase  = np.pi / 2
rf_ref_phase = 0

flip_ex = rf_ex_angle * np.pi / 180
rf_ex, gz, _ = make_sinc_pulse(flip_angle=flip_ex, system=system, duration=n_t_ex/i_raster_time, slice_thickness=slice_thickness,
                               apodization=0.5, time_bw_product=4, phase_offset=rf_ex_phase, return_gz=True)
# rf_ex, gz, _ = make_gauss_pulse(flip_angle=flip_ex, system=system, duration=n_t_ex/i_raster_time, slice_thickness=slice_thickness,
#                                 apodization=0.5, time_bw_product=4, phase_offset=rf_ex_phase)



# def sigpy_n_seq(flip_angle: float, apodization: float = 0, delay: float = 0, duration: float = 0,
#                 freq_offset: float = 0, center_pos: float = 0.5, max_grad: float = 0, max_slew: float = 0,
#                 phase_offset: float = 0, return_gz: bool = True, slice_thickness: float = 0, system: Opts = Opts(),
#                 time_bw_product: float = 4, pulse_cfg: pulse_opts = pulse_opts(),
#                 use: str = str()) -> Union[SimpleNamespace,
#                                            Tuple[
#                                                SimpleNamespace, SimpleNamespace,
#                                                SimpleNamespace]]:
# def make_gauss_pulse(flip_angle: float, apodization: float = 0, bandwidth: float = 0, center_pos: float = 0.5,
#                      delay: float = 0, duration: float = 0, freq_offset: float = 0, max_grad: float = 0,
#                      max_slew: float = 0, phase_offset: float = 0, return_gz: bool = False, slice_thickness: float = 0,
#                      system: Opts = Opts(), time_bw_product: float = 4,
#                      use: str = str()) -> Union[SimpleNamespace,
#                                                 Tuple[SimpleNamespace, SimpleNamespace, SimpleNamespace]]:
#
# def make_sinc_pulse(flip_angle: float, apodization: float = 0, delay: float = 0, duration: float = 0,
#                     freq_offset: float = 0, center_pos: float = 0.5, max_grad: float = 0, max_slew: float = 0,
#                     phase_offset: float = 0, return_gz: bool = False, slice_thickness: float = 0, system: Opts = Opts(),
#                     time_bw_product: float = 4, use: str = str()) -> Union[SimpleNamespace,
#                                                                            Tuple[SimpleNamespace, SimpleNamespace,
#                                                                                  SimpleNamespace]]:

if calc_duration(rf_ex)- t_exwd > 1e-5:   # allow for time check of rf_ex & slice selection flat time
    t_exwd   = calc_duration(rf_ex)
    n_t_exwd = math.ceil(t_exwd * i_raster_time)

flip_ref = rf_flip[0] * np.pi / 180
rf_ref, gz_ref, _ = make_sinc_pulse(flip_angle=flip_ref, system=system, duration=n_t_ref/i_raster_time, slice_thickness=slice_thickness,
                                apodization=0.5, time_bw_product=4, phase_offset=rf_ref_phase, use='refocusing', return_gz=True)
# rf_ref, gz_ref, _ = make_gauss_pulse(flip_angle=flip_ref, system=system, duration=n_t_ref/i_raster_time, slice_thickness=slice_thickness,
#                                      apodization=0.5, time_bw_product=4, phase_offset=rf_ref_phase, use='refocusing')

if calc_duration(rf_ref) - t_refwd > 1e-5:   # allow for time check of rf_ref & slice selection flat time
    t_refwd   = calc_duration(rf_ref)
    n_t_refwd = math.ceil(t_refwd * i_raster_time)

# Save rf pulses - save '.mat'
if save_flag:
    os.chdir('Github/Reconstruction/T2_EPG/rf_pulses')
    parametersMATexc = {'rf_ex':rf_ex}
    nameRFexcsave    = "rf_pulses_exc{}.mat".format(rf_ex_angle)
    parametersMATref = {'rf_ref': rf_ref}
    nameRFrefsave    = "rf_pulses_ref{}.mat".format(rf_ref_angle)
    savemat(nameRFexcsave, parametersMATexc)
    savemat(nameRFrefsave, parametersMATref)
    print("\n\n ----- rf_pulses saved -----  ")

    if rf_pulses:
        print("\n\n quit")
        quit()

#####
#%% --- 4 - Gradients Slice Selection
#####
gs_ex   = make_trapezoid(channel='z', system=system, amplitude=gz.amplitude, flat_time=n_t_exwd/i_raster_time, rise_time=n_dG/i_raster_time)
gs_ref  = make_trapezoid(channel='z', system=system, amplitude=gs_ex.amplitude, flat_time=n_t_refwd/i_raster_time, rise_time=n_dG/i_raster_time)

# Slice Selection spoiling
ags_ex  = gs_ex.area / 2

amplitudeTest = (ags_ex * (1 + fsp_s)) / (n_t_spex/i_raster_time - n_dG/i_raster_time)
# maxGrad*42.54*1e3> amplitudeTest

gs_spr  = make_trapezoid(channel='z', system=system, area=ags_ex * (1 + fsp_s), duration=n_t_sp/i_raster_time, rise_time=n_dG/i_raster_time)
# gs_spr  = make_trapezoid(channel='z', system=system, area=ags_ex * (1 + fsp_s), duration=calc_duration(rf_ex), rise_time=n_dG/i_raster_time)

gs_spex = make_trapezoid(channel='z', system=system, area=ags_ex * fsp_s, duration=n_t_spex/i_raster_time, rise_time=n_dG/i_raster_time)



#####
#%% --- 5 - ADCs / Readouts
#####
# Readout gradient and ADC
gr_acq  = make_trapezoid(channel='x', system=system, flat_area=k_width, flat_time=n_readout_time/i_raster_time, rise_time=n_dG/i_raster_time) # Gradient Readout
n_gr_acq_flat_time = math.ceil(gr_acq.flat_time * i_raster_time)
adc = make_adc(num_samples=Nx, duration=n_gr_acq_flat_time/i_raster_time - 40e-6, delay=20e-6)


#####
#%% --- 6 - Spoilers - Gradient Readout
#####
# RO spoiling
gr_spr  = make_trapezoid(channel='x', system=system, area=gr_acq.area * fsp_r, duration=n_t_sp/i_raster_time, rise_time=n_dG/i_raster_time)
gr_spex = make_trapezoid(channel='x', system=system, area=gr_acq.area * (1 + fsp_r), duration=n_t_spex/i_raster_time, rise_time=n_dG/i_raster_time)

agr_spr = gr_spr.area

# Prephasing gradient in RO direction
agr_preph = gr_acq.area / 2 + agr_spr

#test_Slrate =  ( (agr_preph - k_width)  / (dG**2) ) / gamma  # T/m/s

gr_preph = make_trapezoid(channel='x', system=system, area=agr_preph, duration=n_t_spex/i_raster_time, rise_time=n_dG/i_raster_time)


#####
#%% --- 7 - Defining Phase Areas
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
#%% --- 8 - Split gradients and recombine into blocks - Gradient_Readout & Gradient_Slice_Selection
#####
# gs1 : ramp up of gs_ex
gs1_times = [0, gs_ex.rise_time]
gs1_amp = [0, gs_ex.amplitude]
n_gs1_times = [0, math.ceil(gs_ex.rise_time * i_raster_time)]
n_gs1_amp = [0, math.ceil(gs_ex.amplitude * i_raster_time)]
    
gs1 = make_extended_trapezoid(channel='z', times=[n_gs1_times[0]/i_raster_time, n_gs1_times[1]/i_raster_time], amplitudes=[n_gs1_amp[0]/i_raster_time, n_gs1_amp[1]/i_raster_time])

# gs2 : flat part of gs_ex
gs2_times = [0, gs_ex.flat_time]
gs2_amp = [gs_ex.amplitude, gs_ex.amplitude]
gs2 = make_extended_trapezoid(channel='z', times=gs2_times, amplitudes=gs2_amp)

# gs3 : Bridged slice pre-spoiler
gs3_times = [0, gs_spex.rise_time, gs_spex.rise_time + gs_spex.flat_time,
             gs_spex.rise_time + gs_spex.flat_time + gs_spex.fall_time]
gs3_amp = [gs_ex.amplitude, gs_spex.amplitude, gs_spex.amplitude, gs_ref.amplitude]
gs3 = make_extended_trapezoid(channel='z', times=gs3_times, amplitudes=gs3_amp)

# gs4 : Flat slice selector for pi-pulse
gs4_times = [0, gs_ref.flat_time]
gs4_amp = [gs_ref.amplitude, gs_ref.amplitude]
gs4 = make_extended_trapezoid(channel='z', times=gs4_times, amplitudes=gs4_amp)
calc_duration(gs4)
# gs5 : Bridged slice post-spoiler
gs5_times = [0, gs_spr.rise_time, gs_spr.rise_time + gs_spr.flat_time,
             gs_spr.rise_time + gs_spr.flat_time + gs_spr.fall_time]
gs5_amp = [gs_ref.amplitude, gs_spr.amplitude, gs_spr.amplitude, 0]
gs5 = make_extended_trapezoid(channel='z', times=gs5_times, amplitudes=gs5_amp)
calc_duration(gs5)

# gs7 : The gs3 for next pi-pulse
gs7_times = [0, gs_spr.rise_time, gs_spr.rise_time + gs_spr.flat_time,
             gs_spr.rise_time + gs_spr.flat_time + gs_spr.fall_time]
gs7_amp = [0, gs_spr.amplitude, gs_spr.amplitude, gs_ref.amplitude]
gs7 = make_extended_trapezoid(channel='z', times=gs7_times, amplitudes=gs7_amp)

# gr3 : pre-readout gradient
gr3 = gr_preph

# gr5 : Readout post-spoiler
gr5_times = [0, gr_spr.rise_time, gr_spr.rise_time + gr_spr.flat_time,
             gr_spr.rise_time + gr_spr.flat_time + gr_spr.fall_time]
gr5_amp = [0, gr_spr.amplitude, gr_spr.amplitude, gr_acq.amplitude]
gr5 = make_extended_trapezoid(channel='x', times=gr5_times, amplitudes=gr5_amp)

# gr6 : Flat readout gradient
gr6_times = [0, readout_time]
gr6_amp = [gr_acq.amplitude, gr_acq.amplitude]
gr6 = make_extended_trapezoid(channel='x', times=gr6_times, amplitudes=gr6_amp)

# gr7 : the gr3 for next repeat
gr7_times = [0, gr_spr.rise_time, gr_spr.rise_time + gr_spr.flat_time,
             gr_spr.rise_time + gr_spr.flat_time + gr_spr.fall_time]
gr7_amp = [gr_acq.amplitude, gr_spr.amplitude, gr_spr.amplitude, 0]
gr7 = make_extended_trapezoid(channel='x', times=gr7_times, amplitudes=gr7_amp)

if FBP_test:

    # --- gr_N_2 : N+2 ---
    rf_ex_neg, gz_ex_neg, _ = make_sinc_pulse(flip_angle=-1 * flip_ex, system=system, duration=n_t_ex / i_raster_time,
                                      slice_thickness=slice_thickness,
                                      apodization=0.5, time_bw_product=4, phase_offset=-rf_ex_phase)

    gs_ex_neg  = make_trapezoid(channel='z', system=system, amplitude=gz_ex_neg.amplitude, flat_time=n_t_exwd / i_raster_time,rise_time=n_dG / i_raster_time)
    gs_ref_N_1 = make_trapezoid(channel='z', system=system, amplitude=gs_ex_neg.amplitude, flat_time=n_t_refwd / i_raster_time,rise_time=n_dG / i_raster_time)

    # --- gr_N_1 : N+1 ---
    flip_ref_N_1              = rf_flip_N_1 * np.pi / 180
    rf_ref_N_1, gz_ref_N_1, _ = make_sinc_pulse(flip_angle=flip_ref_N_1, system=system, duration=n_t_ref / i_raster_time,
                                        slice_thickness=slice_thickness,
                                        apodization=0.5, time_bw_product=4, phase_offset=rf_ref_phase, use='refocusing', return_gz=True)

    gs1_end = make_extended_trapezoid(channel='z', times=[n_gs1_times[0] / i_raster_time, n_gs1_times[1] / i_raster_time],
                                  amplitudes=[n_gs1_amp[1] / i_raster_time, n_gs1_amp[0] / i_raster_time])

#####
#%% --- 9 - Calculate timing and delays
#####
# delay_TR : delay at the end of each MSE pulse train (i.e. each TR)
t_ex = gs1.t[-1] + gs2.t[-1] + gs3.t[-1]
t_ref = gs4.t[-1] + gs5.t[-1] + gs7.t[-1] + readout_time
t_end = gs4.t[-1] + gs5.t[-1]
TE_train = t_ex + n_echo * t_ref + t_end
TR_fill = (TR - n_slices * TE_train)
TR_fill = system.grad_raster_time * round(TR_fill / system.grad_raster_time)
if TR_fill < 0:
    TR_fill = 1e-3
    warnings.warn(f'TR too short, adapted to include all slices to: {1000 * n_slices * (TE_train + TR_fill)} ms')
else:
    print(f'TR not filled = {int(1000 * TR_fill)} ms of TR = {int(1000 * TR)} ms, TR fill = {int(1000 * n_slices * TE_train)} ms')


# Delay
delay_TR = make_delay(TR_fill)

# set times in points
n_TE_train = math.ceil(TE_train * i_raster_time)
n_TR_fill  = math.ceil(TR_fill * i_raster_time)
n_delay_TR = math.ceil(delay_TR.delay * i_raster_time)

# interleaved slice
a_list       = list(range(1, n_slices+1))
aux_a        = a_list[::2]
if len(a_list)>1:
    aux_b        = a_list[1::2]
else:
    aux_b = []
n_sli_interl = aux_a + aux_b



#####
#%% --- 10 - Parallel Imaging w/ SENSE factor R
#####
if SENSE:
    n_phA = np.round(int(n_phA/R))
    phase_areas = phase_areas[::R]

if GRAPPA:
    # phase_areas     = pe_steps * delta_k
    # n_phA           = len(phase_areas)
    # restLin         = (Nx-fullLin)/R
    # aux_centerKspac = fullLin/2
    # iniFullcent     = int(n_phA / 2 - aux_centerKspac - 1)                      # Ini of Fully Sampled
    # endFullcent     = int(n_phA / 2 + aux_centerKspac - 1)                      # End of Fully Sampled
    #
    # n_phA           = int(fullLin + ( (Nx - fullLin)/R ))
    # Ny_real         = n_phA                                                     # real Ny for sequence
    #
    # aux_pA_1        = phase_areas[0:iniFullcent:R]                              # Get phases of AF for 1st interval
    # iniKfull        = int(aux_pA_1.shape[0])                                    # Point where it start fully sampled k-space
    # aux_pA_2        = phase_areas[iniFullcent+1:endFullcent+1]                      # Get phases of AF for 2nd interval
    # endKfull        = iniKfull + int(aux_pA_2.shape[0])                         # Point where it ends fully sampled k-space
    # aux_pA_3        = phase_areas[endFullcent+2:int(phase_areas.shape[0]):R]    # Get phases of AF for 3rd interval
    #



    phase_areas     = pe_steps * delta_k
    n_phA           = len(phase_areas)
    restLin         = (Nx-fullLin)/R
    aux_centerKspac = fullLin/2
    iniFullcent     = int(n_phA / 2 - aux_centerKspac - 1)                      # Ini of Fully Sampled
    endFullcent     = int(n_phA / 2 + aux_centerKspac - 1)                      # End of Fully Sampled                                                # real Ny for sequence

    aux_pA_1        = phase_areas[0:iniFullcent:R]                              # Get phases of AF for 1st interval

    aux_val2        = np.where(phase_areas == aux_pA_1[-1])
    iniFullcent_new = int(aux_val2[0]) + 2                                      # Point where it start fully sampled k-space
    endFullcent_new = iniFullcent_new + fullLin                                 # Point where it ends fully sampled k-space
    aux_pA_2        = phase_areas[iniFullcent_new:endFullcent_new]              # Get phases of AF for 2nd interval
    iniKfull        = int(aux_pA_1.shape[0])                                    # Point where it start fully sampled k-space
    endKfull        = iniKfull + int(aux_pA_2.shape[0])                         # Point where it ends fully sampled k-space

    aux_pA_3        = phase_areas[endFullcent_new:int(n_phA):R]     # Get phases of AF for 3rd interval

    n_phA           = int(fullLin + ( (Nx - fullLin)/R ))
    Ny_real         = n_phA

    phase_areas     = np.concatenate((aux_pA_1,aux_pA_2,aux_pA_3),axis=0)

    # n_phA = int(phase_areas.shape[0]) # TESTE 11-4-22 para suprimir o numero de phase areas menor devido ao facto de ter retirado uma linha do kspace


#####
#%% --- 10 - Define sequence blocks/Readouts + Create '.seq' file
#####

# SAR_check = True
#
# while SAR_check:
for ph_A in range(n_phA):  # cicle per TR each with one phase (line of k_space) and multiple Echos
    for s in range(0, len(n_sli_interl)): # For each slice (interleaved)
        rf_ex.freq_offset   = gs_ex.amplitude * slice_thickness * (n_sli_interl[s] - (n_slices - 1) / 2)    # frequency offset for rf pulse - excitation
        rf_ref.freq_offset  = gs_ref.amplitude * slice_thickness * (n_sli_interl[s] - (n_slices - 1) / 2)   # frequency offset for rf pulse - reference
        rf_ex.phase_offset  = rf_ex_phase - 2 * np.pi * rf_ex.freq_offset * calc_rf_center(rf_ex)[0]        # Phase offset for rf pulse - excitation
        rf_ref.phase_offset = rf_ref_phase - 2 * np.pi * rf_ref.freq_offset * calc_rf_center(rf_ref)[0]     # Phase offset for rf pulse - reference

        seq.add_block(gs1)
        seq.add_block(gs2, rf_ex)
        seq.add_block(gs3, gr3)

        for k_echo in range(n_echo): # For each TR
            phase_area = phase_areas[ph_A]
            # Make phase encoding gradients (they are reused for all readouts)
            gp_pre = make_trapezoid(channel='y', system=system, area=phase_area,  duration=n_t_sp/i_raster_time, rise_time=n_dG/i_raster_time)
            gp_rew = make_trapezoid(channel='y', system=system, area=-phase_area, duration=n_t_sp/i_raster_time, rise_time=n_dG/i_raster_time)
            seq.add_block(gs4, rf_ref)
            seq.add_block(gs5, gr5, gp_pre)

            # get ADC
            # if k_ex > 0:
            #     seq.add_block(gr6, adc)
            # else:
            #     seq.add_block(gr6)
            seq.add_block(gr6, adc)

            seq.add_block(gs7, gr7, gp_rew)

        if FBP_test: # Flip Back Pulse - Test
            rf_ex_neg.freq_offset   = gs_ex_neg.amplitude * slice_thickness * (n_sli_interl[s] - (n_slices - 1) / 2)    # frequency offset for rf pulse - excitation
            rf_ref_N_1.freq_offset  = gs_ref_N_1.amplitude * slice_thickness * (n_sli_interl[s] - (n_slices - 1) / 2)   # frequency offset for rf pulse - reference
            rf_ex_neg.phase_offset  = rf_ref_phase - 2 * np.pi * rf_ex_neg.freq_offset * calc_rf_center(rf_ex_neg)[0]        # Phase offset for rf pulse - excitation
            rf_ref_N_1.phase_offset = rf_ref_phase - 2 * np.pi * rf_ref_N_1.freq_offset * calc_rf_center(rf_ref_N_1)[0]     # Phase offset for rf pulse - reference

            # N+1 Pulse --> Teta
            seq.add_block(gs4, rf_ref_N_1)

            # N+2 Pulse --> 90º
            seq.add_block(gs3, gr3)
            seq.add_block(gs2, rf_ex_neg)
            seq.add_block(gs1_end)

        else:
            seq.add_block(gs4)
            seq.add_block(gs5)

    seq.add_block(delay_TR)



# Duration of sequence in ms
seqDur = seq.duration()
print("\n\n .. seq duration .. ", seqDur[0], "in s")


#####
#%% --- 12 - B1+rms - control for not exceeding the allowed radiation for the System
#####
B1plus_rf_ex   = np.max(rf_ex.signal)/gamma * 1e6                     # units (uT)
B1plus_rf_ref  = np.max(rf_ref.signal)/gamma  * 1e6                   # units (uT)
b1plus_t_ex    = flip_ex / (np.max(rf_ex.signal) * 2*np.pi) * 1e3     # time for specific area:  gamma*2pi*B1*time = flip - units (ms)
b1plus_t_refoc = flip_ref / (np.max(rf_ref.signal) * 2*np.pi) * 1e3   # units (ms)

sumb1plus_rms = (  (B1plus_rf_ex ** 2 ) * b1plus_t_ex   +   n_echo * ( (B1plus_rf_ref ** 2) * b1plus_t_refoc )   )  * n_slices * Ny
# sumb1plus_rms = (  (B1plus_rf_ex * b1plus_t_ex)**2      +   n_echo * ( (B1plus_rf_ref*b1plus_t_refoc)**2 )       )  * n_slices * Ny
b1Plus_rms    = np.sqrt( sumb1plus_rms / (seqDur[0]*1e3) )   # units (uT)

print("\n\n .. B1 + rms .. ", b1Plus_rms, "in uT")

if b1Plus_rms>9:
   raise ValueError('B1+rms is higher than threshold - Increase TR')



#####
#%% --- 11 - Plots & Report & K-traject
#####
if plotTest:
    seq.plot()
    # seq.plot(time_range=[0, 0.02 * TR])

if reportTest:
    print(seq.test_report())
    seq.check_timing()

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
#%% --- 12 - save file '.seq'
#####
if save_flag:
    dir = 'Github/Sequence/MSE'
    os.chdir(dir)

    folder_name = 'test%s_MSE'  % (test)
    seq_name = 'test%s_MSE_nEchos-%s_nSlices-%s_rfFlipAng-%s_TE-%sms_TR-%sms_FOV-%smm_Nx-%s_Ny-%s_gm-%s_sm-%s_sliceT-%smm' % (test, n_echo, n_slices, int(rf_flip[0]), round(TE*1e3), round(TR*1e3), round(fov*1e3), Nx, Ny, maxGrad, maxSlew, round(slice_thickness*1e3))
    if os.path.isdir(folder_name):
        print("\n\n - Directory already exist - ")
    else:
        os.makedirs(folder_name)
    os.chdir(folder_name)

    # get definitions
    #seq.set_definition('FOV', [fov, fov, slice_thickness])
    #seq.set_definition('Name', 'MSE_test')

    seq.write(seq_name)

    # Add definitions
    #seq.definitions()

    # save '.mat' file with parameters
    if SENSE: # for Parallel Imaging - SENSE
        parametersMAT = {'test':test, 'TR':TR, 'TE':TE, 'DT':DT, 'nslices':n_slices, 'st':slice_thickness,
                         'max_grad':maxGrad, 'max_slew':maxSlew, 'flipAngle':rf_ex_angle, 'rf_flipAngle':rf_flip,
                         'nEcohs':n_echo, 'Nx':Nx, 'Ny':Ny, 'Ny_real':Ny_real,'FOV':fov, 'duration':seqDur[0],
                         'delta_k':delta_k, 'k_width':k_width, 'R':R}
    elif GRAPPA:  # for Parallel Imaging - SENSE
        parametersMAT = {'test': test, 'TR': TR, 'TE': TE, 'DT': DT, 'nslices': n_slices, 'st': slice_thickness,
                         'max_grad': maxGrad, 'max_slew': maxSlew, 'flipAngle': rf_ex_angle, 'rf_flipAngle': rf_flip,
                         'nEcohs': n_echo, 'Nx': Nx, 'Ny': Ny, 'Ny_real': Ny_real, 'FOV': fov,
                         'duration': seqDur[0], 'delta_k': delta_k, 'k_width': k_width, 'R': R,
                         'fullLin': fullLin, 'iniKfull':iniKfull,'endKfull':endKfull}
    else: # for normal acquisition
        parametersMAT = {'test':test, 'TR':TR, 'TE':TE, 'DT':DT, 'nslices':n_slices, 'st':slice_thickness,
                         'max_grad':maxGrad, 'max_slew':maxSlew, 'flipAngle':rf_ex_angle, 'rf_flipAngle':rf_flip,
                         'nEcohs':n_echo, 'Nx':Nx, 'Ny':Ny, 'FOV':fov, 'duration':seqDur[0], 'delta_k':delta_k,
                         'k_width':k_width}
    savemat("sequence_info.mat", parametersMAT)


    print("\n\n ----- Seq file saved -----  ", seq_name)



#####git
#%% --- 13 - Run seq2xml & Save '.xml'
#####
if ConvertXML:
    os.chdir(dir)

    folder_name = 'test%s_MSE' % (test)
    seq_name = 'test%s_MSE_nEchos-%s_nSlices-%s_rfFlipAng-%s_TE-%sms_TR-%sms_FOV-%smm_Nx-%s_Ny-%s_gm-%s_sm-%s_sliceT-%smm' % (test, n_echo, n_slices, int(rf_flip[0]), round(TE*1e3), round(TR*1e3), round(fov*1e3), Nx, Ny, maxGrad, maxSlew, round(slice_thickness*1e3))
    if os.path.isdir(folder_name):
        print("\n\n - Directory already exist - ")
    else:
        os.makedirs(folder_name)
    os.chdir(folder_name)
    out_folder = dir + '/' + folder_name
    seq2xml(seq, seq_name, out_folder)

    print("\n\n ----- Seq2jemris saved -----  ", seq_name)
