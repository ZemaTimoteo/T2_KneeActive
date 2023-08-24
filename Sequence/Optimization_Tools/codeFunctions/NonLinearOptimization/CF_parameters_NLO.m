function [params] = CF_parameters_NLO(params)

% Other parameters
params.res    = 256;            % Resolution
params.nsli   = 30;
params.B1     = 1;              % B1 value - TODO change
params.T1     = 1000*1e-3;      % T1 (s)
params.T2     = 45*1e-3;        % T2 (s)
params.SNR    = 40;             % TODO - automatizar
params.sigma1 = 2.7200*1e-3;    % RFexcDur in (s) - taken from pulseq implementation of a 90º RF exc pulse
params.sigma2 = 7.7900*1e-3;    % t_gs4 + t_gs5 + t_spoiler in (s) - time of spoiler gradients
params.sigma3 = (1/params.SNR); % sigma of SNR

% directories
params.dict.file_path = 'D:\Tiago\Trabalho\2021_2025_PhD\Projects\qMRI_Joint\Code\matlabCode\qMRI_tools\Sequence_Optimization';
params.dict.dir_data  = [params.dict.file_path '\Data'];
params.dict.dir_rf    = [params.dict.dir_data '\rf_pulses' ];

% Parameters for derivatives - Initialize
params.deriv.R1   = 1/(params.T1);               % T1 in (s)
params.deriv.R2   = 1/(params.T2);               % T2 in (s)
rf_exc            = pi/2;
params.alpha_RF   = rf_exc; % Flip Angle in (rad)

end