function [params] = CF_Mycode_parameters_NLO(params)

% Other parameters
params.res       = 200;                   % Resolution
params.af        = 2;                     % GRAPPA Acceleration Factor
params.accFactor = 1/3 + (2/3)/params.af; % Factor of Acceleration Factor - Full Centre K-space + af of the other 2/3
params.st        = 2.5e-3;                % Slice Thickness (m)
params.nsli      = 25;                    % Number of slices
params.B1        = 1;                     % B1 value - TODO change
params.T1        = 1000;                  % T1 (ms)
params.T2        = 45;                    % T2 (ms)
params.TRini     = 10;                    % Initial TR (ms)
params.SNR       = 40;                    % TODO - automatizar
params.sigma1    = 2.7200;                % RFexcDur in (ms) - taken from pulseq implementation of a 90º RF exc pulse
params.sigma2    = 7.7900;                % t_gs4 + t_gs5 + t_spoiler in (ms) - time of spoiler gradients
params.sigma3    = (1/params.SNR);        % sigma of SNR

% directories & Tests
params.dict.file_path = 'D:\Tiago\Trabalho\2021_2025_PhD\Projects\qMRI_Joint\Code\matlabCode\qMRI_tools\Sequence_Optimization';
params.dict.dir_data  = [params.dict.file_path '\Data'];
params.dict.dir_rf    = [params.dict.dir_data '\rf_pulses' ];
params.plotTest       = 'Fals';

% Parameters for derivatives - Initialize
rf_exc            = pi/2;
params.alpha_RF   = rf_exc; % Flip Angle in (rad)

% Constrains
params.constr.betaMin    = 8;     % Beta min (ms)
params.constr.alphaMin   = 0;     % Alpha min (rad)
params.constr.alphaMax   = pi;    % Alpha max (rad)

end