function [params] = CF_vFA_PS_parameters(params)

% Other parameters
params.B1             = 1;                     % B1 value - TODO change
params.T1             = 1000;                  % T1 (ms)
params.nsli           = 15;                    % Number of slices
params.TotalSlices    = 16;                    % Number of Total of slots for slices (16 due to bit reverse)
params.res            = 230;                   % Resolution
% params.st             = 1.8e-3;                % Slice Thickness (m) - TODO redudant but still here not to crash the code
params.st             = 3e-3;                % Slice Thickness (m) - TODO redudant but still here not to crash the code
% params.t_ex           = 2.93e-3;               % time of exc pulse in (s) - 2.5e-3;
% params.t_ref          = 1.52e-3;               % time of ref pulse in (s) - 1.7e-3;
% params.t_ex           = 2.81e-3;               % time of exc pulse in (s) - 2.5e-3;
% params.t_ref          = 1.6e-3;               % time of ref pulse in (s) - 1.7e-3;
params.t_ex           = 2.8e-3;               % time of exc pulse in (s) - 2.5e-3;
params.t_ref          = 1.45e-3;               % time of ref pulse in (s) - 1.7e-3;
params.npoints        = 25;                    % Number of points for the slr profile
% params.TBWP           = 4;                     % Time Bandwidth Path
params.TBWP_ex        = 4;                     % Time Bandwidth Path
params.TBWP_ref       = 2;                     % Time Bandwidth Path
params.center_pos     = 0.85;                   % Center position of RFexc


params.alpha_RF       = 90*pi/180;
params.TRini          = 10;                     % Initial TR (ms)
params.SNR            = 100;                    % TODO - automatizar
params.sigma1         = params.t_ex;           % RFexcDur in (ms) - taken from pulseq implementation of a 90º RF exc pulse
% t_gs4                 = 0.00162/2*1e3;         % units (ms)
% t_gs5                 = 0.00119*1e3;           % units (ms)
% t_spoiler             = 5;                     % units (ms)
% t_gs4                 = 0.00128/2*1e3;         % units (ms)
% t_gs5                 = 0.00151*1e3;           % units (ms)
% t_gs4                 = 0.00137/2*1e3;         % units (ms)
% t_gs5                 = 0.002*1e3;             % units (ms)
t_gs4                 = 0.00152/2*1e3;         % units (ms)
t_gs5                 = 0.00185*1e3;             % units (ms)
t_spoiler             = 0;                     % units (ms)
params.sigma2         = t_gs4+t_gs5+t_spoiler; % t_gs4 + t_gs5 + t_spoiler in (ms) - time of spoiler gradients
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
    params.af             = 8;   % acceleration factor LORAK
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
params.dict.file_path = params.dir;
params.dict.dir_data  = [params.dict.file_path filesep 'Data'];
params.dict.dir_rf    = [params.dict.dir_data filesep 'rf_pulses' ];
params.plotTest       = 'Fals';

% Parameters for derivatives - Initialize
rf_exc            = pi/2;
params.alpha_exc  = rf_exc; % Flip Angle in (rad)

% Constrains
params.constr.betaMin   = 8;         % TE min (ms)
% params.constr.betaMin   = 6;           % TE min (ms)
params.constr.betaMax   = 13.8;          % TE max (ms)
params.constr.alphaMin  = 1*pi/180;   % FA min (rad)
params.constr.alphaMax  = pi;          % FA max (rad)
params.constr.maxB1_rms = 2.8;           % Max B1rms for SAR (uT)
% params.constr.maxB1_rms = 4.8;           % Max B1rms for SAR (uT)
params.constr.T_acq     = 7;           % Max Time for sequence without acc in min - 10/12min
% params.constr.T_acq     = 9;           % Max Time for sequence without acc in min - 10/12min

% Min Area Under the Curve - from cFA for T2=8ms|T2=45ms|T2=200ms
T2test                  = find([8 45 200 83 30]==params.T2); % To select a minAUC
ETL_cFA                 = 5;
TE_cFA                  = 13.8;

% aux_minAUC              = [1.0161   31.2638  112.9359   61.7597];   % values from cFA
if params.methodDic == 'JUSTdict'  % for dirac
%     aux_minAUC              = [2.70   27.2112 28.2638 20.263];   % values from cFA      4.1359
%     aux_minAUC              = [1.85   28.2 28.26 24.76];   % values from Siemens MSE      4.1359
    % aux_minAUC              = [1.69  21.28 28.26 24.76];   % values from Siemens MSE      4.1359
    aux_minAUC              = [1.69  21.28 28.26 24.76 15];   % values from Siemens MSE      4.1359
%     aux_minAUC              = [ 33.80 21.86*(TE_cFA * ETL_cFA) 1989.4471    785.58323 ];   % AUC/DeltaT | values from cFA      4.1359   
%     aux_maxCRLB           	= [1.13    0.169   0.046  0.071];
    % aux_maxCRLB           	= [0.72  0,0629 0.046  0.071];
    aux_maxCRLB           	= [0.72  0,0629 0.046  0.071 0.004];
elseif params.methodDic == 'SLR_Prof' % For slr profile
    aux_minAUC              = [ 43.97 295.13 1989.44  500.64 ];   % values from cFA      4.1359
    aux_maxCRLB             = [313.20  14.79   70.99   19.09];  % values from cFA      4.1359
end
% aux_minAUC              = [0 0 0 0];   % values from cFA      4.1359

params.constr.minAUC    = aux_minAUC(T2test);
params.constr.minVarT2  = aux_maxCRLB(T2test);      % T2=8ms; ETL=30; TE=40ms - Min Sensitivity of the test varT2 = 0.0027

if T2test == 4 % Test for GM in 3T
    params.T1 = 1445;
end

end