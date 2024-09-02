%% Test CRLB and AUC Values for constant FA optimization
% Test also Time of Scan, SAR, EPG
% by; TTFernandes, Oct - 2023
% @IST

clear all 
clc

%% 0 - Settings
% ... 0.1 - Set Paths ...
addpath(genpath('D:\Tiago\Trabalho\2021_2025_PhD\Projects\qMRI_Joint\Code\matlabCode\qMRI_tools\Sequence_Optimization_vF\MSE_vFA_optimization_vAug24'))

addpath(genpath('D:\Tiago\Trabalho\2021_2025_PhD\Projects\qMRI_Joint\Code\matlabCode\Toolboxes\pulseq-master'));
addpath(genpath('D:\Tiago\Trabalho\2021_2025_PhD\Projects\qMRI_Joint\Code\matlabCode\qMRI_tools\Sequence_Optimization\aux_Functions\rf_tools'));


%% Get parameters
T2_test = 45;
params          = [];
params.Rayleigh      = 'True';                 % Test with noise 'True' or 'Fals'
params.ParallelImg   = 'GRAPP';                % for GRAPPA 'GRAPP' | for LORAKS 'LORAK'
params.RFqual        = 'True';                 % Quality of RF (different slice Thickness in excitation and refocusing - Kumar et al. 2016 JMRI)
testDict             = 'SLR_Prof';             % SLR profile | For EPG & CRLB - 'SLR_Prof' OR 'JUSTdict' OR '90DirSLR' OR '90SLRDir'
params.methodDic     = testDict;               % For EPG & CRLB - 'SLR_Prof' OR 'JUSTdict' OR '90DirSLR' OR '90SLRDir'
params.T2            = T2_test;                % T2 (ms) [8 45 200 83]
params.dir_rf        = 'D:\Tiago\Trabalho\2021_2025_PhD\Projects\qMRI_Joint\Code\matlabCode\qMRI_tools\Sequence_Optimization_vF\MSE_vFA_optimization_vAug24\Data\rf_pulses';

% params          = CF_Mycode_parameters_NLO(params); % Set parameters for Cost Function (CF)
params          = CF_PS_parameters(params); % Set parameters for Cost Function (CF)

display(['Dictionary tested with ' testDict])

%% Test cFA CRLB and AUC - from ISMRM 23 abstract TTFernandes - CRLB minimization for optimal Multi Spin-Echo T2-mapping of the cartilage using dictionary-based estimation
% parameters
% % T1_ini  = [ones(1,3)*1000 1445];
% % T2_ini  = [   8   45  200   83];
% % TE_ini  = [   8    8 8.14 8.03];
% % ETL_ini = [  10   16   26   26];
% % FA_ini  = [ 151  148  147  146];
% % TR_ini  = [2387 3528 5622 5549];

% GRAPPA 
% % T1_ini  = [ones(1,3)*1000 1445];
% % T2_ini  = [   8   45  200   83];
% % TE_ini  = [8.25 8.25 8.14    9];
% % ETL_ini = [  10   12   26   12];
% % FA_ini  = [ 153  152  147  158];
% % TR_ini  = [2394 2807 5622 3041];

% last collumn is theoretic from mean of references
T1_ini  = [ones(1,2)*1000 1445  1000];
T2_ini  = [   8   45        83    45];
TE_ini  = [8.25 8.25         9 10.17];
ETL_ini = [  10   12        12     7];
FA_ini  = [ 153  152       158   180];
TR_ini  = [2394 2807      3041  2309];

% % % LORAKS 
% % T1_ini  = [ones(1,3)*1000 1445];
% % T2_ini  = [   8   45  200   83];
% % TE_ini  = [8.25 8.25 8.14  8.3];
% % ETL_ini = [  10   12   26   26];
% % FA_ini  = [ 153  150  147  149];
% % TR_ini  = [4086 3632 5622 4898];

display(['T2 values used for cFA optimiz. (ms)           : ' num2str(T2_ini)])

%% test CRLB
CRLBini = zeros(1,size(T1_ini,2));
for jj=1:size(T1_ini,2)
    params.T1    = T1_ini(jj);
    params.T2    = T2_ini(jj);
    params.TRini = TR_ini(jj);
    x_ini{jj}    = [TE_ini(jj) FA_ini(jj)*ones(1,ETL_ini(jj))*pi/180];      

    CRLBini(jj) = Constraint_test_CRLB(x_ini{jj},ETL_ini(jj),params);    % check sensitivity to estimate T2 values
end
display(['Values for CRLB for cFA optimiz.               : ' num2str(CRLBini)])

%% testAUC
AUCini  = zeros(1,size(T1_ini,2));
for jj=1:size(T1_ini,2)
    params.T1    = T1_ini(jj);
    params.T2    = T2_ini(jj);
    params.TRini = TR_ini(jj);

    AUCini(jj) = Constraint_test_AUC(x_ini{jj},ETL_ini(jj),params);    % check sensitivity to estimate T2 values
end
display(['Values for AUC for cFA optimiz.                : ' num2str(AUCini)])

%% epg figure
figure()
for jj=1:size(T1_ini,2)
    params.T1    = T1_ini(jj);
    params.T2    = T2_ini(jj);
    params.TRini = TR_ini(jj);
    
    data       = CF_epg_optFSE_getSignalEPG(x_ini{jj},ETL_ini(jj),params);
    signal{jj} = data.signal;
        
    plot(signal{jj},'linewidth',2)
    hold on
end
title(['EPG of each T2 optz. for constant Flip Angles'])
% % legend(['Short | 8    (ms)'],['Knee | 45   (ms)'],['GM | 200 (ms)'],['GM   | 83   (ms)'])
legend(['Short | 8    (ms)'],['Knee | 45   (ms)'],['GM   | 83   (ms)'],['Knee (theory) | 45   (ms)'])

%% Get Time of Scan
TRini  = zeros(1,size(T1_ini,2));
ToSini  = zeros(1,size(T1_ini,2));
for jj=1:size(T1_ini,2)
    params.T1    = T1_ini(jj);
    params.T2    = T2_ini(jj);
    params.TRini = TR_ini(jj);

    [TRini(jj), ToSini(jj)] = Constraint_test_Time(x_ini{jj},ETL_ini(jj),params);    % check sensitivity to estimate T2 values
end
display(['Values for Time of Scan for cFA optimiz. (min) : ' num2str(ToSini)])
display(['Values for TR for cFA optimiz. (s)             : ' num2str(TRini)])

%% Get SAR
SARini  = zeros(1,size(T1_ini,2));
for jj=1:size(T1_ini,2)
    params.T1    = T1_ini(jj);
    params.T2    = T2_ini(jj);
    params.TRini = TR_ini(jj);

    SARini(jj) = Constraint_test_SAR(x_ini{jj},ETL_ini(jj),params);    % check sensitivity to estimate T2 values
end
display(['Values for SAR(uT) for cFA optimiz.            : ' num2str(SARini)])

