function [vardT2] = CF_CRLB_epg_optFSE(x,params)
% Obtain the variance expected for a specific model/set of parameters
%
% Functions used:
%   slr_profile
%   epg_derivatives_vF: to generate the states over the slice excited (sum) for different T2
%
% Inputs:
%   x = [TE (ms), FA(rad)]
%
% Ouputs:
%   test_vardT2: Variance with respect to T2 estimation.
%

%% ... 1 - Parameters ...

% 1.1 - Inputs
TE    = x(1);              % beta = TE - Echo Time
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

% 1.3 - directories & tests
plotTest   = 'Fals';
methodDic  = 'JUSTdict';

% 1.4 - Parameters for derivatives - Initialize
beta     = TE*1e-3;         % interecho spacing in (s)
tau      = beta;        % Full of echospacing (ESP) in (ms) = TE
FAngles  = FA;              % Flip Angle in (rad)
alpha_RF = params.alpha_RF; % Flip Angle in (rad)


%% 2 - Get theta and d_theta
% 2.1 - theta formulation
numer_theta_beta = 1 - exp( - ( ( sigma1/2 + (1/2+ETL)*beta )*(nsli-1)+sigma2*nsli) / T1 );
denom_theta_beta = (sigma3^2/res) * sqrt( (sigma1/2 + ( 1/2 + ETL )*beta + sigma2)*nsli*res);
theta.theta_beta = (numer_theta_beta) ^ 2 / denom_theta_beta;


%% ... 3 - Get Gradient ...
% 3.1 - Calculate derivatives - % (from toolbox T2estimation-master) and adpted

[ ~, grad, ~ ] = epg_derivatives_vF( ETL, tau, R1, R2, alpha_RF, FAngles);

ds_dT2 = grad.ds_dT2; % derivative over T2

%% ... 4 - CRLB ...
    % --- 4.1 Uncertainty of CRLB ---
vardT2 = theta.theta_beta * (ds_dT2(:)'*ds_dT2(:));

end

