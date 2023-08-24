function [vardT2, grad] = CF_Mycode_Grad_CRLB_epg_optFSE_NLO(x,params)
% Cost Function and Gradient obtained for the variance expected for a specific model/set of parameters
%
% Functions used:
%   CF_Mycode_epg_derivatives_NLO: to generate the states over the slice excited (sum) for different T2
%
% Inputs:
%   x = [TE (ms), FA(rad)]
%   params - parameters
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
T1     = params.T1;      % T1 (ms)
T2     = params.T2;      % T2 (ms)
SNR    = params.SNR;     
sigma1 = params.sigma1;  % RFexcDur in (ms) - taken from pulseq implementation of a 90º RF exc pulse
sigma2 = params.sigma2;  % t_gs4 + t_gs5 + t_spoiler in (ms) - time of spoiler gradients
sigma3 = params.sigma3;  % sigma of SNR

% 1.3 - directories
plotTest   = params.plotTest;
file_path  = 'D:\Tiago\Trabalho\2021_2025_PhD\Projects\qMRI_Joint\Code\matlabCode\qMRI_tools\Sequence_Optimization';
dir_data   = [file_path '\Data'];

% 1.4 - Parameters for derivatives - Initialize
beta     = TE;              % interecho spacing in (ms) - Full of echospacing (ESP) in (ms) = TE
alpha_RF = params.alpha_RF; % Flip Angle in (rad)


%% 1 - Get gamma and d_gamma
% 1.1 - gamma thetformulation
numer_gamma_beta = (  1 - exp( - ( ( sigma1/2 + (1/2+ETL)*beta )*(nsli-1)+sigma2*nsli) / T1 ) )  ^ 2;
denom_gamma_beta = (sigma3^2/res) * sqrt( (sigma1/2 + ( 1/2 + ETL )*beta + sigma2)*nsli*res);
gamma.gamma_beta = numer_gamma_beta / denom_gamma_beta;

% 1.2 - Derivative of gamma in order to beta (ms)
d_numer_gamma_beta_dbeta = (  2*(1/2 + ETL)    *    (-1 + nsli)  *  exp(  (  -nsli*sigma2 - (-1 + nsli)*( sigma1/2 + (1/2 + ETL)*beta ) ) / T1  )  * ...
                              ( 1 - exp(  (-nsli*sigma2 - (-1 + nsli)*(sigma1/2 + (1/2 + ETL)*beta))/T1  ) ) ...
                             ) /  T1;
d_denom_gamma_beta_dbeta = ((1/2 + ETL) * nsli * sigma3^2) / (  sqrt(2)*sqrt(nsli*res * (sigma1 + 2*sigma2 + beta + 2*ETL*beta))  );


gamma.dgamma_beta = ( denom_gamma_beta * d_numer_gamma_beta_dbeta  - numer_gamma_beta * d_denom_gamma_beta_dbeta) / ...
                         (denom_gamma_beta^2)  ;

                     
%% ... 2 - Calculate Gradients (derivatives) of EPG
aux_grad = CF_Mycode_epg_derivatives_NLO(ETL, beta, ...
    T1, T2, alpha_RF, FA, params, gamma); % [dF_dalphaL,dF_dbeta]

ds_dT2 = aux_grad.ds_dT2; % derivative over T2


%% ... 3 - CRLB ...
    % --- 3.1 Calculate CF value - Uncertainty of CRLB ---
vardT2 = - gamma.gamma_beta * (ds_dT2(:)'*ds_dT2(:)); % multiplied by (-) in order to maximize
% % vardT2 = gamma.gamma_beta * (ds_dT2(:)'*ds_dT2(:)); 

%% ... 4 - Gradients
if nargout > 1
    grad = [aux_grad.dF_dbeta aux_grad.dF_dalpha];   % Gradient of [dF_dbeta, dF_dalphaL] 
end

end

