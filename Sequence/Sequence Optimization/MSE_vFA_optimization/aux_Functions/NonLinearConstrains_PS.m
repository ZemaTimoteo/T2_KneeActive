%% File details
%
% Nonlinear Constraints for optimization problem of Fast Spin Echo (FSE) T2
% Search Pattern
%   Mapping. Take into account 3 constraints:
%          - SAR
%          - Time of Scan
%          - AUC
%
% by: TTFernandes, Octorber 2023
% @IST, Uni Lisbon

function [c,ceq] = NonLinearConstrains_PS(x,ETL,params)

addpath(genpath('D:\Tiago\Trabalho\2021_2025_PhD\Projects\qMRI_Joint\Code\matlabCode\Toolboxes\pulseq-master'));

%% 1 - Import a seq file to compute SAR for
% ... 1.05 - Parameters from input ...
st_exc    = params.st_ex;          % Slice Thickness (m) - '2.6e-3'
st_ref    = params.st_ref;          % Slice Thickness (m) - '2.6e-3'
accFactor = params.accFactor;   % Factor of Acceleration for GRAPPA or LORAKS
nsli      = params.nsli;        % number of slices - '30'
res       = params.res;         % Resolution Nx = Ny
new_res   = res*accFactor;      % Resolution Needed
TR        = params.TRini;       % (ms)

TE	 = x(1);    %(ms)
FA   = x(2:end);

% ... 1.1 - Parameters for pulseq ...
maxGrad        = 32;                % value for 1.5T & 3T - gm =  30| 32 mT/m,  respectively
maxSlew        = 120;               % value for 1.5T & 3T - sm = 124|130 T/m/s, respectively
rf_ex_phase    = pi/2;              % RF excitation Phase
rf_ref_phase   = 0;                 % RF refocusing Phase
flip_ex        = params.alpha_RF;   % Flip angle excitation in rad
t_ex           = params.t_ex;       % in ms
t_ref          = params.t_ref;      % in ms 
tend           = 0.0079*1e3;        % Spoilers gradients times in (s)
gamma          = 42.54e6;           % Gyromagnetic constant (Hz/T)


% ... 1.2 - initialize sequence ...
seq = mr.Sequence();
sys = mr.opts('MaxGrad', maxGrad, 'GradUnit', 'mT/m', 'MaxSlew', maxSlew, ...
              'SlewUnit', 'T/m/s', 'rfRingdownTime', 100e-6, ...
              'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6);

% ... 1.3 - Partial derivatives ....
deriv_B1plus_rf_ref_dFApartial  = zeros(1,ETL);
deriv_b1plus_t_refoc_dFApartial = zeros(1,ETL);
deriv_sumb1plus_rms_dFApartial  = zeros(1,ETL);
deriv_b1Plus_rms_dFApartial     = zeros(1,ETL);


%% 2 - Identify RF blocks and compute SAR - 10 seconds must be less than twice and 6 minutes must be less than 4 (WB) and 3.2 (head-20)
T_vector   = [];
SARwbg_vec = zeros(1,size(FA,2));
SARhg_vec  = zeros(1,size(FA,2));

% ... 2.1 - Impact of  rf_excitation ...
% Create rf_excitation
[rf_exc, gz]   = mr.makeSincPulse(flip_ex,'Duration',t_ex,...
    'SliceThickness',st_exc,...
    'apodization',0.5,'timeBwProduct',4,...
    'phaseOffset',rf_ex_phase,'system',sys);
rf_excDur      = mr.calcDuration(rf_exc);
T_vector       = rf_excDur/2*1e3;
signal         = rf_exc.signal;

    % Calculate maxSignal amplitude - Global
B1plus_rf_ex   = (max(signal)/gamma * 1e6);                 % (uT)
b1plus_t_ex    = flip_ex /  (max(signal) * 2*pi) * 1e3;    % time for specific area (1/uT): trf = flip_angle / b1_max


% ... 2.2 - Impact of rf_refocusing ...
for jj=1:size(FA,2)
    % Create rf_refocusing
    flip_ref       = FA(jj);  % Flip angle refocusing in rad
    [rf_refoc, gz] = mr.makeSincPulse(flip_ref, 'Duration', t_ref,...
                                    'SliceThickness',st_ref,...
                                    'apodization',0.5, 'timeBwProduct',4,...
                                    'phaseOffset',rf_ref_phase, 'use','refocusing',...
                                    'system',sys);
    rf_refocDur = mr.calcDuration(rf_refoc);

    % Time vector
    if jj == 1
        T_vector = [T_vector T_vector(jj)+(TE/2)];
    else
        T_vector = [T_vector T_vector(jj)+TE];
    end
    
    signal    = rf_refoc.signal;
    
    % % flipDeriv       = rf_refoc.flipDeriv;
    % % aux_signalDeriv = rf_refoc.signalDeriv;
    % % signalDeriv     = aux_signalDeriv/flipDeriv;
    
    % Calculate maxSignal amplitude - Global
    B1plus_rf_ref(jj)            = max(signal)/gamma * 1e6;                                       % (uT)
    b1plus_t_refoc(jj)           = (FA(jj)) / (max(signal) * 2*pi) * 1e3;    % time for specific area (1/uT): trf = flip_angle / b1_max

    % Calculate derivatives for each FA partial
    deriv_B1plus_rf_ref_dFApartial(jj)  = 0;                                                     % because of the maximum function
    deriv_b1plus_t_refoc_dFApartial(jj) = (max(signal) + 0) / ((max(signal))^2 * 2*pi) * 1e3;    % time for specific area (s): trf = flip_angle / b1_max    
end

T_vector = [T_vector T_vector(end)+tend];


%% 3 - Obtain Time vector
RFexcDur   = rf_excDur *1e3;    % (ms)
RFrefocDur = rf_refocDur *1e3;  % (ms)

if (TE/2- RFexcDur/2 - RFrefocDur/2)<0 || (TE - RFrefocDur)<0
    b1Plus_rms  = NaN;
    T_scan_s    = NaN;
    return
end

% crushers
t_gs4     = 0.00198/2*1e3; % units (ms)
t_gs5     = 0.0018*1e3;    % units (ms)
t_spoiler = 5;             % units (ms)

% Get time
aux_Trec  = RFexcDur + (TE/2- RFexcDur/2) + TE * ETL;  % (ms)
aux_TRmin = aux_Trec + t_gs4 + t_gs5 + t_spoiler;      % (ms)
TRmin     = aux_TRmin * nsli;                          % (ms)  - Wait for all slices to be filled out

% Respect T2 map image - recovery of longitudinal magnt big enought.
if TRmin<TR 
    T_scan  = TR * new_res;                   % (ms)
    Trec    = TR - aux_Trec;                        % (ms)
    TR      = TR*1e-3;                              % (s)
else
    T_scan  = TRmin * new_res;                         % (ms) - Time of Scan
    Trec    = TRmin - aux_Trec;                     % (ms) - Recuperation Time
    TR      = TRmin*1e-3;                           % (s)  - Repetition Time
end

T_scan_s   = T_scan*1e-3;                              % (s)
T_scan_m   = T_scan_s/60;                              % (min)


%% 4 - Time averaged RF power - match Siemens data
% ... 4.1 - Average RFpower
sumb1plus_rms = ( (B1plus_rf_ex )^2 * b1plus_t_ex +  sum(B1plus_rf_ref.^2 .*b1plus_t_refoc ) ) * ...
                    nsli * new_res;
                
b1Plus_rms    = sqrt(sumb1plus_rms/T_scan);    % units (uT) - eq.2+1/2 fo Keerthivasan, 2019


%% 5 - Get constrains
% 5.1 - Constrains
% Max B1_rms (uT) | maxB1_rms = 5;   % uT
maxB1_rms = params.constr.maxB1_rms;   % uT

% Max Time for GRAPPA % LORAKS (s) | T_acq     = 12;                 % Time in min - 10/12min
maxTime   = params.constr.T_acq;       % Full k-space centre + 2/3 accelerated by af - units(min)
maxTime_s = maxTime * 60;              % Time in (s)
    
% 5.2 - They need to be minor than previous settled constrains
% Inequal constrains - NEED to BE ALL NEGATIVE 
c(1) = b1Plus_rms - maxB1_rms;              % check if it is under the limits for SAR
c(2) = T_scan_s - maxTime_s;                % check if it is under a specific time of scan
c(3) = Constraint_PS_AUC(x,ETL,params);   % check constraintAUC - needs to be bigger than standart - but still negative

% Equal constrains (not needed for this problem)
ceq  = [];


