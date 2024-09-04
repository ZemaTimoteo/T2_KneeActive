function [params] = CF_CRLB_parameters(B1,T1,nslices,res,sliceThickn,TR,SNR,params)

% Other parameters
params.B1             = B1;                    % B1 value - TODO change
params.T1             = T1;                    % T1 (ms)
params.nsli           = nslices;               % Number of slices
params.res            = res;                   % Resolution
params.st             = sliceThickn.exc;       % Slice Thickness (m) - TODO redudant but still here not to crash the code
params.t_ex           = 2.81e-3;               % time of exc pulse in (s) - 2.5e-3;
params.t_ref          = 1.6e-3;                % time of ref pulse in (s) - 1.7e-3;

params.center_pos     = 0.5;                   % Center position of RFexc

params.alpha_exc      = pi/2;                   % Excitation Flip Angle in (rad)
params.TRini          = TR;                     % Initial TR (ms)
params.SNR            = SNR;                    

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


end