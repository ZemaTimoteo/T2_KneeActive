function [vardT2,grad] = CRLB_epg_optFSE_fmincon_TRvar(x)
% Obtain the variance expected for a specific model/set of parameters
%
% Functions used:
%   slr_profile
%   epg_derivatives_vF: to generate the states over the slice excited (sum) for different T2
%
% Inputs:
%   x = [TE (ms), TR (ms), ETL, FA(º)]
%
% Ouputs:
%   CRLB_num:    CRLB value
%   uCRLB:       CRLB uncertainty (Std_CRLB / s) - in Zhang et al. 2013 in Proc. EMBS
%   test_vardT2: Variance with respect to T2 estimation.
%


%% ... 0 - Parameters ...

% Inputs
dTE   = x(1);               % Echo Time
% % ETL   = x(2);               % Echo Train Lenght
% % FA    = x(3);               % Flip Angle
ETL    = 6;                % Echo Train Lenght
FA     = x(2)*ones(1,ETL);  % Flip Angle constant
TRacq  = 3.1162;            % TR of acquisition (s)
Trec   = 2.99939e+03;       % Time of Recovery (ms) - TODO - Valor vem da funçao de constrains

% Other parameters
res    = 256;          % Resolution
T_scan = TRacq*res;    % Time of acquisition (s)
rf_exc = pi/2;
B1     = 1;            % B1 value - TODO change
T1     = 1000;         % T1 (ms)
T2     = 45;           % T2 (ms)
SNR    = 40;           % TODO - automatizar
sigma  = (1/SNR);


plotTest   = 'Fals';
testFBP    = 'Fals';
file_path  = 'D:\Tiago\Trabalho\2021_2025_PhD\Projects\qMRI_Joint\Code\matlabCode\qMRI_tools\Sequence_Optimization';
dir_data   = [file_path '\Data'];
methodDic  = 'JUSTdict';


%% ... 1 - Derivative (from toolbox T2estimation-master) ...
% 1.1 - Initialize
TE     = dTE*1e-3;                  % interecho spacing in (s)
R1     = 1/(T1*1e-3);               % T1 in (s)
R2     = 1/(T2*1e-3);               % T2 in (s)
dir_rf = [dir_data '\rf_pulses' ];

% 1.2 - Calculate derivatives
if methodDic == 'SLR_Prof'    
    clear grad_tool FF_tool
    % 1.2.1 - RF pulses
    refoc_puls = [];
    for jj=1:size(FA,2)
        [exc_puls, aux_refoc_puls] = slr_profile_vF(rf_exc,B1,FA(jj),dir_rf);
        refoc_puls                 = [refoc_puls aux_refoc_puls(:)];  vv
    end

    % 1.2.2 - Get Derivatives
    [ ~, grad_tool, FF_tool ] = epg_derivatives_fmincon_vF(ETL,TE,R1,R2,exc_puls, refoc_puls);
    
elseif methodDic == 'JUSTdict' 
    clear grad_tool FF_tool
    % 1.2.1 - RF pulses    
    FA_exc     = 90;
    FAngles    = FA*pi/180;   % Flip Angle in (rad)
    alpha_RF   = FA_exc*pi/180;  % Flip Angle in (rad)
    
    % 1.2.2 - Get Derivatives
    [ ~, grad_tool, FF_tool ] = epg_derivatives_fmincon_vF(ETL,TE,R1,R2,alpha_RF, FAngles);

elseif methodDic == '90DirSLR'
    clear grad_tool FF_tool
    % 1.2.1 - RF pulses
    refoc_puls = [];
    for jj=1:size(FA,2)
        [~, aux_refoc_puls] = slr_profile(B1,FA(jj),dTE,'HSM2',dir_data);
        refoc_puls                 = [refoc_puls aux_refoc_puls(:)];  
    end 
    alpha_RF  = FA_exc*pi/180;  % Flip Angle in (rad)            
    exc_pulse = repmat(alpha_RF,size(refoc_pulse,1),1); 
    
    % 1.2.2 - Get Derivatives
    [ ~, grad_tool, FF_tool ] = epg_derivatives_fmincon_vF(ETL,TE,R1,R2,exc_pulse, refoc_puls);

elseif methodDic == '90SLRDir'
    clear grad_tool FF_tool
    % 3.2.1 - RF pulses
    refoc_puls = [];
    exc_puls   = [];
    for jj=1:size(FA,2)
        [exc_puls, ~] = slr_profile(B1,FA(jj),dTE,'HSM2',dir_data);
    end    
    radAngle   = FA*pi/180*B1;  % Flip Angle in (rad)
    refoc_puls = repmat(radAngle,size(exc_pulse_3,1),size(FA,2));
    
    % 1.2.2 - Get Derivatives
    [ ~, grad_tool, FF_tool ] = epg_derivatives_fmincon_vF(ETL,TE,R1,R2,exc_pulse,refoc_puls);  
    
end

% 1.3 - Get Variables
ds_dT2 = grad_tool.ds_dT2; % derivative over T2
ds_dFA = grad_tool.ds_dFA; % derivative over Flip Angle


   % --- Matrix derivatives for fmincon ---
% 1.35 - gradients from epg_derivatives_fmincon
ds_dT2_dTE  = grad_tool.ds_dT2_dTE;
ds_dT2_dETL = grad_tool.ds_dT2_dETL;
ds_dT2_dFA  = grad_tool.ds_dT2_dFA;
   % --------------------------------------


% 1.4 - normalize ??? fará sentido?
% % norm_epgData     = FF_tool ./ norm(FF_tool);
% % norm_ds_tool_dT2 = ds_dT2  ./ norm(ds_dT2);

norm_epgData     = FF_tool ;
norm_ds_tool_dT2 = ds_dT2;

% 1.5 ========= Control for FBP - (eq.3 - Keerthivasan, 2019) =============
if testFBP == 'True'
    % 1.5.1 - Control for Trec - EPG
    aux_norm_epgData = norm_epgData;
    auxVariab_S_T2 = exp( - Trec/T1); % control for Trecovery & TR variation
    aux2var_S_T2   = (1 - auxVariab_S_T2) / ( 1 - auxVariab_S_T2 * norm_epgData(ETL) );
    aux2var_S_dT2  = (1 - auxVariab_S_T2) / ( 1 - auxVariab_S_T2 * norm_ds_tool_dT2(ETL) );
    norm_epgData   = aux2var_S_T2.* norm_epgData;   

    % 1.5.2 - Calculate factor
    factorFPB_epg = mean(abs(norm_epgData./aux_norm_epgData));
    
    % 1.5.3 - Control for Trec - Derivative
    norm_ds_tool_dT2 = aux2var_S_dT2.*norm_epgData + aux2var_S_T2.* norm_ds_tool_dT2;   
else
    factorFPB_epg    = 1 - exp( - Trec /T1); % control for Trecovery & TR variation
    norm_epgData     = factorFPB_epg.* norm_epgData;   
    norm_ds_tool_dT2 = factorFPB_epg.* norm_ds_tool_dT2;   

    
       % --- Matrix derivatives for fmincon ---    
    norm_ds_tool_dT2_dTE  = factorFPB_epg.* ds_dT2_dTE;   
    norm_ds_tool_dT2_dETL = factorFPB_epg.* ds_dT2_dETL;   
    norm_ds_tool_dT2_dFA  = factorFPB_epg.* ds_dT2_dFA;   
       % --------------------------------------    

end


%% ... 2 - Figures ...
if plotTest == 'True'   
    % 2.1 - Plot EPGs
    figure()
    plot(abs(norm_epgData),'r'), hold on  % Toolbox EPG
    legend('EPGs Toolb States')
    title(['Method: ',methodDic,', RFexc = ',num2str(rf_exc),', RFref = ',num2str(FA(1))])
    
    % 2.2 - Plot dEPGs
    figure()
    plot(abs(norm_ds_tool_dT2), 'gx--'), hold on
    legend('dT2 EPG Toolbox')            
end

%% ... 3 - Jacobian Matrix ...
% % dM_num  = [dS_dM0(:) norm_ds_tool_dT2(:)]; 
dS_dT2_dEPG = norm_ds_tool_dT2;

   % --- Matrix derivatives for fmincon ---    
dS_dT2_dTE  = norm_ds_tool_dT2_dTE;   
dS_dT2_dETL = norm_ds_tool_dT2_dETL;   
dS_dT2_dFA  = norm_ds_tool_dT2_dFA;   
   % --------------------------------------    
       
%% ... 4 - CRLB ...

	% --- 4.1 CRLBnum ---
% % FIM_num    = dM_num.'   *  (   (1/sigma^2) * eye( size(dM_num,1) )  )   *  dM_num; % Fisher matrix
% % FIM_numT2  = FIM_num(2,2);
% % CRLB_num   = diag(  ( inv(FIM_num) )  );

    % --- 4.2 Uncertainty of CRLB ---
aux_vardT2 = ( (dS_dT2_dEPG(:)'*dS_dT2_dEPG(:)) / (sigma^2) ) / (TRacq); % Divided by TR? 

    % --- 4.3 Time balanced uCRLB by scan time ---
% It is minus 1, because the goal is to maximize this value    
vardT2 = (-1) * aux_vardT2*sqrt(T_scan);   % Divided by the sqrt of TimeScan(s)

%% ... 5 - Gradients ...
if nargout > 1 % gradient required
    aux1 =  sqrt(T_scan) / ( (sigma^2) * (TRacq) );
    grad_vardT2_dTE  = (  dS_dT2_dTE(:)'*dS_dT2_dEPG(:)  +  dS_dT2_dEPG(:)'*dS_dT2_dTE(:)  ) * aux1; % d_vardT2 / d_dTE
%     grad_vardT2_dETL = (  dS_dT2_dETL(:)'*dS_dT2_dEPG(:) +  dS_dT2_dEPG(:)'*dS_dT2_dETL(:) ) / aux1; % d_vardT2 / d_dETL
    grad_vardT2_dFA  = (  dS_dT2_dFA(:)'*dS_dT2_dEPG(:)  +  dS_dT2_dEPG(:)'*dS_dT2_dFA(:)  ) * aux1; % d_vardT2 / d_dFA
    
    % Gradients for fmincon
%     grad = [grad_vardT2_dTE; grad_vardT2_dETL; grad_vardT2_dFA];
    grad = [grad_vardT2_dTE; grad_vardT2_dFA];
end

end

