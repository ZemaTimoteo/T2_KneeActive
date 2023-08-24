function [vardT2, grad] = CF_Grad_CRLB_epg_optFSE_NLO(x,params)
% Cost Function and Gradient obtained for the variance expected for a specific model/set of parameters
%
% Functions used:
%   slr_profile
%   epg_derivatives_vF: to generate the states over the slice excited (sum) for different T2
%
% Inputs:
%   x = [TE (ms), FA(rad)]
%
% Ouputs:
%   vardT2: Variance with respect to T2 estimation.
%   grad: Gradients of each variable according to NLO annotations



%% ... 1 - Parameters ...

% 1.1 - Inputs
TE    = x(1);              % beta = TE - Echo Time (ms)
ETL   = params.ETL;        % K = Echo Train Lenght = ETL
FA    = x(2:end);          % Flip Angle constant (rad) 

% 1.2 - Other parameters
res    = params.res;     % Resolution
nsli   = params.nsli;    % Number of slices
B1     = params.B1;      % B1 value - TODO change
T1     = params.T1;      % T1 (s)
T2     = params.T2;      % T2 (s)
SNR    = params.SNR;     
sigma1 = params.sigma1;  % RFexcDur in (s) - taken from pulseq implementation of a 90º RF exc pulse
sigma2 = params.sigma2;  % t_gs4 + t_gs5 + t_spoiler in (s) - time of spoiler gradients
sigma3 = params.sigma3;  % sigma of SNR

% 1.3 - directories
plotTest   = 'Fals';
file_path  = 'D:\Tiago\Trabalho\2021_2025_PhD\Projects\qMRI_Joint\Code\matlabCode\qMRI_tools\Sequence_Optimization';
dir_data   = [file_path '\Data'];
dir_rf     = [dir_data '\rf_pulses' ];
methodDic  = 'JUSTdict';

% 1.4 - Parameters for derivatives - Initialize
beta     = TE*1e-3;              % interecho spacing in (ms)
tau      = beta;            % Full of echospacing (ESP) in (ms) = TE
R1       = params.deriv.R1; % T1 in (ms^-1)
R2       = params.deriv.R2; % T2 in (ms^-1)
FAngles  = FA;              % Flip Angle in (rad)
alpha_RF = params.alpha_RF; % Flip Angle in (rad)


%% 1 - Get theta and d_theta
% 1.1 - theta formulation
numer_theta_beta = (  1 - exp( - ( ( sigma1/2 + (1/2+ETL)*beta )*(nsli-1)+sigma2*nsli) / T1 ) )  ^ 2;
denom_theta_beta = (sigma3^2/res) * sqrt( (sigma1/2 + ( 1/2 + ETL )*beta + sigma2)*nsli*res);
theta.theta_beta = numer_theta_beta / denom_theta_beta;

% 1.2 - Derivative of theta in order to beta (ms)
d_numer_theta_beta_dbeta = (  2*(1/2 + ETL)    *    (-1 + nsli)*exp(  (  -nsli*sigma2 - (-1 + nsli)*( sigma1/2 + (1/2 + ETL)*beta ) ) / T1  )  * ...
                              ( 1 - exp(  (-nsli*sigma2 - (-1 + nsli)*(sigma1/2 + (1/2 + ETL)*beta))/T1  ) ) ...
                             ) /  T1;
d_denom_theta_beta_dbeta = ((1/2 + ETL) * nsli * sigma3^2) / (  sqrt(2)*sqrt(nsli*res * (sigma1 + 2*sigma2 + beta + 2*ETL*beta))  );


theta.dtheta_beta = ( denom_theta_beta * d_numer_theta_beta_dbeta  - numer_theta_beta * d_denom_theta_beta_dbeta) / ...
                         (denom_theta_beta^2)  ;


%% ... 2 - Get Cost Function and Gradient ...
% (from toolbox T2estimation-master) and adpted
% 2.1 - Calculate derivatives
 grad = CF_epg_derivatives_NLO(ETL, tau, ...
               R1, R2, alpha_RF, FAngles, theta); % [dF_dalphaL,dF_dbeta]

 ds_dT2 = grad.ds_dT2; % derivative over T2

 
%% Teste
% % % 2.1 - Calculate derivatives
% %  grad = CF_Mycode_epg_derivatives_NLO(ETL, tau, ...
% %                R1, R2, alpha_RF, FAngles, theta); % [dF_dalphaL,dF_dbeta]
% % 
% %  ds_dT2 = grad.ds_dT2; % derivative over T2

 
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
% % fprintf(['vardT2 = ' num2str(vardT2) '\n'])

end

