% Tests SE sequence for T2 estimation for mono-exponential protocol dicom data for phantom (HLuz)
% TTFernandes 16_1_2024
% @IST

%% 0 - Set matlab paths and toolboxes
clear all
clc

addpath(genpath('D:\Tiago\Trabalho\2021_2025_PhD\Projects\qMRI_Joint\Code\matlabCode\MSE_recon'));
addpath(genpath('D:\Tiago\Trabalho\2021_2025_PhD\Projects\GIRFOS_tool\Code\0_matlabCode\10_toolboxes'))
addpath(genpath('D:\Tiago\Trabalho\2021_2025_PhD\Projects\qMRI_Joint\Code\matlabCode\Toolboxes\Read_dat'))
%% 0.5 - Tests
tic

plotTest    = 'True';
saveResults = 'True';
maskTest    = 'True';       % 'True' - Get mask             OR 'Fals' - load mask
readPython  = 'Fals';       % 'getIt' - Get Dictionary      OR 'ready' - load Dictionary
% tempMatch = 'getIt';      % 'getIt' - Get ind_parameters  OR 'ready' - load ind_parameters

sliceTest = [1];
numSli    = size(sliceTest,2);   % Number of Slices

toc
%% 1 - load Data SE

% ..... 1.1 - get dir ...
% myCD = 'Z:\Datasets\T2KneeActive_project';
myCD = 'C:\Users\filia\Documents\PhD\Projetos\qMRI\Sequences\MSE\24_01_23_Sequences\test344_MSE_cFA_T2_45ms_test2';
cd(myCD)

str                  = {'select file'};
[file_chos,data_dir] = uigetfile;
cd(data_dir)
twix_obj   = mapVBVD(file_chos);
twix_obj{1,2}.hdr.Protocol.tFree

