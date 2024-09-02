function [params] = CF_PS_parameters(params)

% Other parameters
params.B1             = 1;                     % B1 value - TODO change
params.T1             = 1000;                  % T1 (ms)
params.nsli           = 25;                    % Number of slices
params.res            = 256;                   % Resolution
params.st             = 2.5e-3;                % Slice Thickness (m) - TODO redudant but still here not to crash the code
% params.t_ex           = 2.93e-3;               % time of exc pulse in (s) - 2.5e-3;
% params.t_ref          = 1.52e-3;               % time of ref pulse in (s) - 1.7e-3;
params.t_ex           = 2.81e-3;               % time of exc pulse in (s) - 2.5e-3;
params.t_ref          = 1.6e-3;               % time of ref pulse in (s) - 1.7e-3;
params.npoints        = 52;                    % Number of points for the slr profile

params.center_pos     = 0.5;                   % Center position of RFexc

params.alpha_RF       = 90*pi/180;
params.TRini          = 10;                    % Initial TR (ms)
params.SNR            = 40;                    % TODO - automatizar
params.sigma1         = params.t_ex;           % RFexcDur in (ms) - taken from pulseq implementation of a 90º RF exc pulse
params.sigma2         = 7.7900;                % t_gs4 + t_gs5 + t_spoiler in (ms) - time of spoiler gradients
params.sigma3         = (1/params.SNR);        % sigma of SNR

% Quality fo RF (different slice Thickness in excitation and refocusing - Kumar et al. 2016 JMRI)
if params.RFqual == 'True'
    params.st_ex          = params.st;             % Slice Thickness (m)
    params.st_ref         = params.st_ex*3;        % Slice Thickness (m)
else
    params.st_ex          = params.st;             % Slice Thickness (m)
    params.st_ref         = params.st;             % Slice Thickness (m)    
end


% Acceleration
% GRAPPA
if params.ParallelImg == 'GRAPP'
    params.af        = 2;        % acceleration factor GRAPPA
    params.accFactor = 1/3 + (2/3)/params.af; % Factor of Acceleration Factor - Full Centre K-space + af of the other 2/3
% LORAKS
elseif params.ParallelImg == 'LORAK' 
    params.af             = 4;   % acceleration factor LORAK
    params.partialFourier = 1;
    params.accFactor      = 1/5 + params.partialFourier * (2/5)/params.af + (2/5)/params.af; % Factor of Acceleration Factor - Full Centre K-space + af of the other 2/3
end

% Noise - Rician Distribution
if params.Rayleigh == 'True'
%     r          = makedist('Rician','s',0,'sigma',params.sigma3); % sigma3 = 1/SNR
%     noise_r      = random(r,1,1);
    r            = makedist('Rayleigh','B',params.sigma3); % sigma3 = 1/SNR
    meanval_r    = params.sigma3*sqrt(pi/2);
    params.noise = meanval_r;    
else
    params.noise = 0;
end

% directories & Tests
params.dict.file_path = 'D:\Tiago\Trabalho\2021_2025_PhD\Projects\qMRI_Joint\Code\matlabCode\qMRI_tools\Sequence_Optimization';
params.dict.dir_data  = [params.dict.file_path '\Data'];
params.dict.dir_rf    = [params.dict.dir_data '\rf_pulses' ];
params.plotTest       = 'Fals';

% Parameters for derivatives - Initialize
rf_exc            = pi/2;
params.alpha_exc  = rf_exc; % Flip Angle in (rad)

% Constrains
params.constr.betaMin   = 8.25;        % TE min (ms)
params.constr.betaMax   = 30;          % TE max (ms)
params.constr.alphaMin  = 60*pi/180;   % FA min (rad)
params.constr.alphaMax  = pi;          % FA max (rad)
params.constr.maxB1_rms = 5;           % Max B1rms for SAR (uT)
params.constr.T_acq     = 9;           % Max Time for sequence without acc in min - 10/12min

% Min Area Under the Curve - from cFA for T2=8ms|T2=45ms|T2=200ms
T2test                  = find([8 45 200 83]==params.T2); % To select a minAUC
% aux_minAUC              = [1.0161   31.2638  112.9359   61.7597];   % values from cFA  
if params.methodDic == 'JUSTdict'  % for dirac
    aux_minAUC              = [1.0161   20.2638 28.2638 20.263];   % values from cFA      4.1359
    aux_maxCRLB           	= [1.1947  	 0.1836	 0.0464  0.123];
elseif params.methodDic == 'SLR_Prof' % For slr profile
    aux_minAUC          = [119.59343   1224.1082  1989.4471    785.58323];   % values from cFA      4.1359
    aux_maxCRLB         = [1773.3092  	241.3822  	70.9949  	52.283074];  % values from cFA      4.1359
end
% aux_minAUC              = [0 0 0 0];   % values from cFA      4.1359

params.constr.minAUC    = aux_minAUC(T2test);
params.constr.minVarT2  = aux_maxCRLB(4);      % T2=8ms; ETL=30; TE=40ms - Min Sensitivity of the test varT2 = 0.0027

if T2test == 4 % Test for GM in 3T
    params.T1 = 1445;
end

end