function [vardT2] = CRLB_epg_optFSE_NLO(x,params)
% Cost Function obtained for the variance expected for a specific model/set of parameters
%
% Functions used:
%   slr_profile
%   epg_derivatives_vF: to generate the states over the slice excited (sum) for different T2
%
% Inputs:
%   x = [TE (ms), TR (ms), ETL, FA(º)]
%
% Ouputs:
%   vardT2: Variance with respect to T2 estimation.
%   grad: Gradients of each variable according to NLO annotations



%% ... 0 - Parameters ...

% Inputs
TE    = x(1);              % beta = TE - Echo Time
ETL   = params.ETL;        % K = Echo Train Lenght = ETL
FA    = x(2:end);  % Flip Angle constant

% Other parameters
res    = 256;          % Resolution
nsli   = 30;
B1     = 1;            % B1 value - TODO change
T1     = 1000;         % T1 (ms)
T2     = 45;           % T2 (ms)
SNR    = 40;           % TODO - automatizar
sigma1 = 2.7200*1e-3;  % RFexcDur in (s) - taken from pulseq implementation of a 90º RF exc pulse
sigma2 = 7.7900*1e-3;  % t_gs4 + t_gs5 + t_spoiler in (s) - time of spoiler gradients
sigma3 = (1/SNR);      % sigma of SNR

% directories
plotTest   = 'Fals';
file_path  = 'D:\Tiago\Trabalho\2021_2025_PhD\Projects\qMRI_Joint\Code\matlabCode\qMRI_tools\Sequence_Optimization';
dir_data   = [file_path '\Data'];
dir_rf     = [dir_data '\rf_pulses' ];
methodDic  = 'JUSTdict';

% Parameters for derivatives - Initialize
beta     = TE*1e-3;                  % interecho spacing in (s)
R1       = 1/(T1*1e-3);               % T1 in (s)
R2       = 1/(T2*1e-3);               % T2 in (s)
rf_exc   = pi/2;
FAngles  = FA*pi/180;     % Flip Angle in (rad)
alpha_RF = rf_exc; % Flip Angle in (rad)
tau      = beta;          % Full of echospacing (ESP) in (s) = TE


%% 1 - Get theta and d_theta
% 1.1 -convert all in seconds
T1_s   = T1*1e-3;

% 1.3 - theta formulation
numer_theta_beta = 1 - exp( - ( ( sigma1/2 + (1/2+ETL)*beta )*(nsli-1)+sigma2*nsli) / T1_s );
denom_theta_beta = (sigma3^2/res) * sqrt( (sigma1/2 + ( 1/2 + ETL )*beta + sigma2)*nsli*res);
theta.theta_beta = (numer_theta_beta) ^ 2 / denom_theta_beta;

% 1.4 - Derivative of theta in order to beta
theta.dtheta_beta = - ( res^2 * (sigma1/2 + (ETL + 1/2)*beta + sigma2) * (1 - exp((-(nsli - 1)*(sigma1/2 + (ETL + 1/2)*beta) - sigma2*nsli) / T1_s))) / ...
                      (  2 * sigma3^2 * (   res*nsli*( sigma1/2 + (ETL + 1/2)*beta + sigma2)  )^(3/2) ) ...
                    - ( res * (-sigma1/2 - (ETL + 1/2)*beta - sigma2) * exp( (-(nsli - 1)*(sigma1/2 + (ETL + 1/2)*beta) - sigma2*nsli) / T1_s )  )  / ...
                      ( sigma3^2 * T1_s * sqrt( res*nsli * (sigma1/2 + (ETL + 1/2)*beta + sigma2)));
            
            
%% ... 2 - Get Cost Function and Gradient ...
% (from toolbox T2estimation-master) and adpted
% 2.1 - Calculate derivatives
 [ obj, grad, FF ] = epg_derivatives_NLO_vF(ETL, tau, ...
     R1, R2, alpha_RF, FAngles, theta); % [dF_dalphaL,dF_dbeta]

 ds_dT2 = grad.ds_dT2; % derivative over T2

%% ... 4 - Figures ...
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

%% ... 5 - CRLB ...
    % --- 4.1 Uncertainty of CRLB ---
vardT2 = theta.theta_beta * (ds_dT2(:)'*ds_dT2(:));
fprintf(['vardT2 = ' num2str(vardT2) '\n'])

end

