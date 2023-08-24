#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon October 25 12:33:25 2021
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
from scipy.io import savemat

PC = 0  # Seia=1 or myPC=0

if PC == 1:
    sys.path.insert(1, '/home/tfernandes/Documents/PYTHON/Toolboxes/pypulseq')
    os.chdir('/home/tfernandes/Documents/PYTHON/Toolboxes/pypulseq')
    from Sequence.sequence import Sequence
    from calc_duration import calc_duration
    from make_adc import make_adc
    from make_delay import make_delay
    from make_sinc_pulse import make_sinc_pulse
    from make_sinc_pulse_channel import make_sinc_pulse_channel
    from make_trap_pulse import make_trapezoid
    from make_arbitrary_grad import make_arbitrary_grad
    from make_extended_trapezoid import make_extended_trapezoid
    from opts import Opts

    sys.path.insert(1, '/home/tfernandes/Documents/PYTHON/Toolboxes/pypulseq/Sequence')

    sys.path.append('/home/tfernandes/Documents/PYTHON/Toolboxes')
    os.chdir('/home/tfernandes/Documents/Projetos/Project_lfStroke/Code/1_pythonCodes/B0_ECC')

elif PC == 0:
    import pypulseq

    from pypulseq.Sequence.sequence import Sequence
    from pypulseq.calc_rf_center import calc_rf_center
    from pypulseq.make_adc import make_adc
    from pypulseq.make_delay import make_delay
    from pypulseq.make_extended_trapezoid import make_extended_trapezoid
    from pypulseq.make_sinc_pulse import make_sinc_pulse
    from pypulseq.make_trap_pulse import make_trapezoid
    from pypulseq.opts import Opts
    from pypulseq.calc_duration import calc_duration
    from pypulseq.make_block_pulse import make_block_pulse

    sys.path.append('D:/Tiago/Trabalho/2021_2025_PhD/Projects/Toolboxes/Python/py2jemris')
    # os.chdir('D:/Tiago/Trabalho/2021_2025_PhD/Projects/Toolboxes/Python/py2jemris')

from seq2xml import seq2xml

# %% ===========================================================================
#####
# %% --- 0.2 - Settings
#####
test = 253

save_flag  = True  # save file
plotTest   = False  # Plot sequence
reportTest = False  # Report the sequence
ConvertXML = False  # Convert seq2xml

#####
# --- 1 - Set sequence limits ---
#####

# scan parameters
maxGrad = 30  # in mT/m
maxSlew = 124  # in T/m/s
dG      = 200e-6  # rise time (s)

n_echo   = 20  # #echoes
n_slices = 1  # #slices

flip_ang = 90  # RF excitation - in degrees
rf_flip = 170  # RF refocusing - in degrees

TE      = 7e-3                   # in s (min: 8e-3)
TR_fill = 0.001                  # in (s)  delay at the end of each MSE pulse train (i.e. each TR)
TR      = (n_echo+1)*TE+TR_fill  # in s
TE_eff  = 60e-3                  # Effective TE in s
pe_type = 'linear'

i_raster_time = 100000
DT = 1 / i_raster_time  # in s

# image parameters
fov = 125e-3  # in m
Nx, Ny = 50, 50
slice_thickness = 2.6e-3  # in m - (min: 2.6e-3)

# aux_calc
k0 = round(TE_eff / TE)
delta_k = 1 / fov
k_width = Nx * delta_k

gamma   = 42.54e6  # Hz/Tesla



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
readout_time = 2.5e-3 + 2 * system.adc_dead_time  # time of readout
t_ex = 2e-3  # time of excitation rf_pulse   - 2.9376e-3
t_exwd = t_ex + system.rf_ringdown_time + system.rf_dead_time  # time of excitation window| time of excitation + ringdown time + dead time
t_ref = 1e-3  # time of refocusing rf_pulse  - 1.5232e-3
t_refwd = t_ref + system.rf_ringdown_time + system.rf_dead_time  # time of refocusing window| time of refocusing + ringdown time + dead time
t_sp = 0.5 * (TE - readout_time - t_refwd)  # time of spoiler| TE>dG & cannot be smaller than readout_time + refwindow
t_spex = 0.5 * (TE - t_exwd - t_refwd)  # time of spoiler excitation|

if TE < 4 * dG + readout_time + t_refwd:
    raise ValueError('TE is too small compared with rise time (dG)')

# timing in number of points
n_TE = math.ceil(TE * i_raster_time)
n_readout_time = math.ceil(readout_time * i_raster_time)
n_t_ex = math.ceil(t_ex * i_raster_time)
n_t_exwd = math.ceil(t_exwd * i_raster_time)
n_t_ref = math.ceil(t_ref * i_raster_time)
n_t_refwd = math.ceil(t_refwd * i_raster_time)
n_t_sp = math.ceil(t_sp * i_raster_time)
n_t_spex = math.ceil(t_spex * i_raster_time)
n_dG = math.ceil(dG * i_raster_time)

fsp_r = 1
fsp_s = 0.5



#####
# %% --- 3 - RF pulses
#####
rf_ex_phase = np.pi / 2
rf_ref_phase = 0

flip_ex = flip_ang * np.pi / 180
rf_ex, gz, _ = make_sinc_pulse(flip_angle=flip_ex, system=system, duration=n_t_ex / i_raster_time,
                               slice_thickness=slice_thickness,
                               apodization=0.5, time_bw_product=4, phase_offset=rf_ex_phase)

if calc_duration(rf_ex) - t_exwd > 1e-5:  # allow for time check of rf_ex & slice selection flat time
    t_exwd = calc_duration(rf_ex)
    n_t_exwd = math.ceil(t_exwd * i_raster_time)

flip_ref = rf_flip[0] * np.pi / 180
rf_ref, gz, _ = make_sinc_pulse(flip_angle=flip_ref, system=system, duration=n_t_ref / i_raster_time,
                                slice_thickness=slice_thickness,
                                apodization=0.5, time_bw_product=2, phase_offset=rf_ref_phase, use='refocusing')

if calc_duration(rf_ref) - t_refwd > 1e-5:  # allow for time check of rf_ref & slice selection flat time
    t_refwd = calc_duration(rf_ref)
    n_t_refwd = math.ceil(t_refwd * i_raster_time)



#####
# %% --- 4 - Gradients Slice Selection
#####
gs_ex = make_trapezoid(channel='z', system=system, amplitude=gz.amplitude, flat_time=n_t_exwd / i_raster_time,
                       rise_time=n_dG / i_raster_time)
gs_ref = make_trapezoid(channel='z', system=system, amplitude=gs_ex.amplitude, flat_time=n_t_refwd / i_raster_time,
                        rise_time=n_dG / i_raster_time)

# Slice Selection spoiling
ags_ex = gs_ex.area / 2

amplitudeTest = (ags_ex * (1 + fsp_s)) / (n_t_spex / i_raster_time - n_dG / i_raster_time)
# maxGrad*42.54*1e3> amplitudeTest

gs_spr = make_trapezoid(channel='z', system=system, area=ags_ex * (1 + fsp_s), duration=n_t_sp / i_raster_time,
                        rise_time=n_dG / i_raster_time)
# gs_spr  = make_trapezoid(channel='z', system=system, area=ags_ex * (1 + fsp_s), duration=calc_duration(rf_ex), rise_time=n_dG/i_raster_time)

gs_spex = make_trapezoid(channel='z', system=system, area=ags_ex * fsp_s, duration=n_t_spex / i_raster_time,
                         rise_time=n_dG / i_raster_time)



#####
# %% --- 5 - ADCs / Readouts
#####
# Readout gradient and ADC
gr_acq = make_trapezoid(channel='x', system=system, flat_area=k_width, flat_time=n_readout_time / i_raster_time,
                        rise_time=n_dG / i_raster_time)  # Gradient Readout
n_gr_acq_flat_time = math.ceil(gr_acq.flat_time * i_raster_time)
adc = make_adc(num_samples=Nx, duration=n_gr_acq_flat_time / i_raster_time - 40e-6, delay=20e-6)



#####
# %% --- 6 - Spoilers - Gradient Readout
#####
# RO spoiling
gr_spr = make_trapezoid(channel='x', system=system, area=gr_acq.area * fsp_r, duration=n_t_sp / i_raster_time,
                        rise_time=n_dG / i_raster_time)
gr_spex = make_trapezoid(channel='x', system=system, area=gr_acq.area * (1 + fsp_r), duration=n_t_spex / i_raster_time,
                         rise_time=n_dG / i_raster_time)

agr_spr = gr_spr.area

# Prephasing gradient in RO direction
agr_preph = gr_acq.area / 2 + agr_spr
gr_preph = make_trapezoid(channel='x', system=system, area=agr_preph, duration=n_t_spex / i_raster_time,
                          rise_time=n_dG / i_raster_time)

#####
# %% --- 7 - Defining Phase Areas
#####
# Phase encoding lobe calculations
pe_steps = np.floor(np.arange(1, Ny + 1) - 0.5 * Ny - 1)
phase_areas = pe_steps * delta_k
n_phA = len(phase_areas)

#####
# %% --- 8 - Split gradients and recombine into blocks - Gradient_Readout & Gradient_Slice_Selection
#####
# gs1 : ramp up of gs_ex
gs1_times = [0, gs_ex.rise_time]
gs1_amp = [0, gs_ex.amplitude]
n_gs1_times = [0, math.ceil(gs_ex.rise_time * i_raster_time)]
n_gs1_amp = [0, math.ceil(gs_ex.amplitude * i_raster_time)]

gs1 = make_extended_trapezoid(channel='z', times=[n_gs1_times[0] / i_raster_time, n_gs1_times[1] / i_raster_time],
                              amplitudes=[n_gs1_amp[0] / i_raster_time, n_gs1_amp[1] / i_raster_time])

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



#####
# %% --- 9 - Calculate timing and delays
#####
# Delay
delay_pEx = make_delay(gs1.t[-1])
delay_hTE = make_delay(TE/2-gs4.t[-1]/2-gs2.t[-1]/2)
delay_TE  = make_delay(TE-gs4.t[-1])
delay_TR  = make_delay(TR_fill)

# set times in points
n_delay_pEx = math.ceil(delay_pEx.delay * i_raster_time)
n_delay_hTE = math.ceil(delay_hTE.delay * i_raster_time)
ñ_delay_TE  = math.ceil(delay_TE.delay * i_raster_time)
n_delay_TR  = math.ceil(delay_TR.delay * i_raster_time)

# interleaved slice
a_list = list(range(1, n_slices + 1))
aux_a = a_list[::2]
if len(a_list) > 1:
    aux_b = a_list[1::2]
else:
    aux_b = []
n_sli_interl = aux_a + aux_b



#####
# %% --- 10 - Define sequence blocks/Readouts + Create '.seq' file
#####
for ph_A in range(n_phA):  # cicle per TR each with one phase and multiple Echos
    for s in range(0, len(n_sli_interl)):  # For each slice (interleaved)
        # print(n_sli_interl[s])
        rf_ex.freq_offset = gs_ex.amplitude * slice_thickness * (
                    n_sli_interl[s] - (n_slices - 1) / 2)  # frequency offset for rf pulse - excitation
        rf_ref.freq_offset = gs_ref.amplitude * slice_thickness * (
                    n_sli_interl[s] - (n_slices - 1) / 2)  # frequency offset for rf pulse - reference
        rf_ex.phase_offset = rf_ex_phase - 2 * np.pi * rf_ex.freq_offset * calc_rf_center(rf_ex)[
            0]  # Phase offset for rf pulse - excitation
        rf_ref.phase_offset = rf_ref_phase - 2 * np.pi * rf_ref.freq_offset * calc_rf_center(rf_ref)[
            0]  # Phase offset for rf pulse - reference

        seq.add_block(delay_pEx)
        seq.add_block(rf_ex)
        seq.add_block(delay_hTE)

        for k_echo in range(n_echo):  # For each TR
            seq.add_block(rf_ref)
            seq.add_block(delay_TE)

    seq.add_block(delay_TR)

# Duration of sequence in ms
seqDur = seq.duration()
print("\n\n .. seq duration .. ", seqDur[0], "in s")



#####
# %% --- 12 - B1+rms - control for not exceeding the allowed radiation for the System
#####
B1plus_rf_ex   = np.max(rf_ex.signal) / gamma * 1e6  # units (uT)
B1plus_rf_ref  = np.max(rf_ref.signal) / gamma * 1e6  # units (uT)
b1plus_t_ex    = flip_ex / (np.max(rf_ex.signal) * 2 * np.pi) * 1e3  # time for specific area:  gamma*2pi*B1*time = flip - units (ms)
b1plus_t_refoc = flip_ref / (np.max(rf_ref.signal) * 2 * np.pi) * 1e3  # units (ms)

sumb1plus_rms = ((B1plus_rf_ex ** 2 ) * b1plus_t_ex+ n_echo * ((B1plus_rf_ref ** 2) * b1plus_t_refoc)) * n_slices * Ny
b1Plus_rms = np.sqrt(sumb1plus_rms / (seqDur[0] * 1e3))  # units (uT)

print("\n\n .. B1 + rms .. ", b1Plus_rms, "in uT")

if b1Plus_rms > 9:
    raise ValueError('B1+rms is higher than threshold - Increase TR')



#####
# %% --- 11 - Plots & Report
#####
if plotTest:
    seq.plot()
    # seq.plot(time_range=[0, 0.02 * TR])

if reportTest:
    print(seq.test_report())
    seq.check_timing()

#####
# %% --- 12 - save file '.seq'
#####
if save_flag:
    if PC == 1:
        os.chdir('/home/tfernandes/Documents/Projetos/Project_Cartilage/Data/qMRI/ETs')
    else:
        dir = 'C:/Users/filia/Documentss/PhD/Projetos/qMRI/Sequences/ETs'
        os.chdir('C:/Users/filia/Documents/PhD/Projetos/qMRI/Sequences/ETs')

    folder_name = 'test%s_ETs' % (test)
    seq_name = 'test%s_ETs_Echos-%s_nSlices-%s_rfFlipAng-%s_TE-%sms_TR-%sms_FOV-%smm_Nx-%s_Ny-%s_gm-%s_sm-%s_sliceT-%smm' % (
    test, n_echo, n_slices, int(rf_flip[0]), round(TE * 1e3), round(TR * 1e3), round(fov * 1e3), Nx, Ny, maxGrad,
    maxSlew, round(slice_thickness * 1e3))
    if os.path.isdir(folder_name):
        print("\n\n - Directory already exist - ")
    else:
        os.makedirs(folder_name)
    os.chdir(folder_name)
    seq.write(seq_name)

    # save '.mat' file with parameters
    parametersMAT = {'test': test, 'TR': TR, 'TE': TE, 'DT': DT, 'nslices': n_slices, 'st': slice_thickness,
                     'max_grad': maxGrad, 'max_slew': maxSlew, 'flipAngle': flip_ang, 'rf_flipAngle': rf_flip,
                     'nEcohs': n_echo, 'Nx': Nx, 'Ny': Ny, 'FOV': fov, 'duration': seqDur[0], 'delta_k': delta_k,
                     'k_width': k_width}
    savemat("sequence_info.mat", parametersMAT)

    print("\n\n ----- Seq file saved -----  ", seq_name)
