

function [vardT2, ds_dT2, FF_tool, factorFPB_epg] = CRLB_epg_optFSE_testAUC(rf_exc, B1, T1, T2, dTE, ETL, TRacq, T_scan, Trec, flipA, sigma, plotTest,testFBP, dir_data)
% Obtain the variance expected for a specific model/set of parameters
%
% Functions used:
%   slr_profile
%   my_epg: to generate the states over the slice excited (sum) for different T2
%
% Inputs:
%   FA_exc:      Flip angle excitation pulse: vector
%   B1_scale:    B1 field
%   T1:          constant- relaxation time of the longitudinal magnetization
%   T2:          interval- relaxation time of the transversal magnetization
%   dTE:         Spacing Time between Echoes (ms)
%   ETL:         refocusing echoes train length
%   Trec:        Recovery time (ms)
%   flipA:       flip angle for refocusing pulse (degrees)
%   sigma:       Noise variance for CRLB
%   step:        step for numerical derivative
%   plotTest:    'True' or 'Fals'
%
% Ouputs:
%   CRLB_num:    CRLB value
%   uCRLB:       CRLB uncertainty (Std_CRLB / s) - in Zhang et al. 2013 in Proc. EMBS
%   test_vardT2: Variance with respect to T2 estimation.
%


%% ... 1 - Derivative (from toolbox T2estimation-master) ...
% 1.1 - Initialize
nEchoes    = ETL;
TE         = dTE*1e-3;      % interecho spacing in (s)
R1         = 1/(T1*1e-3);   % T1 in (s)
R2         = 1/(T2*1e-3);   % T2 in (s)
dir_rf     = [dir_data '\rf_pulses' ];


% 1.2 - Calculate derivatives
clear grad_tool FF_tool
% 1.2.1 - RF pulses
FA_exc     = 90;
flipAngles = flipA*pi/180;   % Flip Angle in (rad)
alpha_RF   = FA_exc*pi/180;  % Flip Angle in (rad)

% 1.2.2 - Get Derivatives and EPG
[ ~, grad_tool, FF_tool ] = epg_derivatives_vF(nEchoes,TE,R1,R2,alpha_RF, flipAngles);


% 1.3 - Get Variables
ds_dT2           = grad_tool.ds_dT2; % derivative over T2
ds_dFA           = grad_tool.ds_dFA; % derivative over Flip Angle
norm_epgData     = FF_tool ;
norm_ds_tool_dT2 = ds_dT2;


% 1.5 ========= Control for FBP - (eq.3 - Keerthivasan, 2019) =============
factorFPB_epg    = 1- exp( - Trec /T1); % control for Trecovery & TR variation - need for specific T1
norm_epgData     = factorFPB_epg.* norm_epgData;
norm_ds_tool_dT2 = factorFPB_epg.* norm_ds_tool_dT2;

%% ... 2 - Get AUC ...
x = (1:ETL)*TE;   % Square brackets waste time here only
A = trapz(x, abs(norm_epgData));

%% ... 3 - Figures ...
if plotTest == 'True'   
    % 2.1 - Plot EPGs
    figure()
    plot(abs(norm_epgData),'r'), hold on  % Toolbox EPG
    legend('EPGs Toolb States')
    title(['Method: ',methodDic,', RFexc = ',num2str(rf_exc),', RFref = ',num2str(flipA(1))])
    
    % 2.2 - Plot dEPGs
    figure()
    plot(abs(norm_ds_tool_dT2), 'gx--'), hold on
    legend('dT2 EPG Toolbox')            
end

%% ... 4 - Jacobian Matrix ...
% % dM_num  = [dS_dM0(:) norm_ds_tool_dT2(:)]; 
dS_dT2_dEPG = norm_ds_tool_dT2;

%% ... 5 - CRLB ...

	% --- 4.1 CRLBnum ---
% % FIM_num    = dM_num.'   *  (   (1/sigma^2) * eye( size(dM_num,1) )  )   *  dM_num; % Fisher matrix
% % FIM_numT2  = FIM_num(2,2);
% % CRLB_num   = diag(  ( inv(FIM_num) )  );

    % --- 4.2 Uncertainty of CRLB ---
% % uCRLB       = (abs(CRLB_num(2))/(T2^2))*(TRacq);    % variance CRLB / s
% - in Zhang et al. 2013 in Proc. EMBS % corrected also for the SNR due to distribution of the noise
aux_vardT2 = ( (dS_dT2_dEPG(:)'*dS_dT2_dEPG(:)) / (sigma^2) ) / (TRacq); % divided by TR in order to balanced the available signal depending on the signal available 

% --- 4.3 Time balanced uCRLB by scan time ---
vardT2 = aux_vardT2*sqrt(T_scan);   % Divided by the sqrt of TimeScan(s)

end

