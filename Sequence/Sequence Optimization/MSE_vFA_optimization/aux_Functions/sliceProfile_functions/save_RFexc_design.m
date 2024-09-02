function save_RFexc_design(params,st_ex,rf_exc,dir_data,name_rf_exc_pulse)

%% 1 - Import a seq file to compute SAR for
% ... 1.1 - Parameters for pulseq ...
maxGrad        = 30;                % value for 1.5T & 3T - gm =  30| 32 mT/m,  respectively
maxSlew        = 120;               % value for 1.5T & 3T - sm = 124|130 T/m/s, respectively
sliceThickness = st_ex;             % slice thickness
rf_ex_phase    = pi/2;              % RF excitation Phase
rf_ref_phase   = 0;                 % RF refocusing Phase
flip_ex        = rf_exc;            % Flip angle excitation in rad
t_ex           = params.t_ex;       % in ms


% ... 1.2 - initialize sequence ...
seq = mr.Sequence();
sys = mr.opts('MaxGrad', maxGrad, 'GradUnit', 'mT/m', 'MaxSlew', maxSlew, ...
              'SlewUnit', 'T/m/s', 'rfRingdownTime', 100e-6, ...
              'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6);


%% 2 - Identify RF blocks and compute SAR - 10 seconds must be less than twice and 6 minutes must be less than 4 (WB) and 3.2 (head-20)
T_vector   = [];

% ... 2.1 - Impact of  rf_excitation ...
% Create rf_excitation
[rf_exc, gz]   = mr.makeSincPulse(flip_ex,'Duration',t_ex,...
    'SliceThickness',sliceThickness,...
    'apodization',0.5,'timeBwProduct',4,...
    'phaseOffset',rf_ex_phase,'system',sys);
rf_excDur      = mr.calcDuration(rf_exc);
T_vector       = rf_excDur/2*1e3;
signal         = rf_exc.signal;

%% 3 - Save RF
% 3.1 - Save RFexcitation
rf_ex.type          = 'rf';
rf_ex.signal        = signal;
rf_ex.t             = rf_excDur;
rf_ex.freq_offset   = 0;
rf_ex.phase_offset  = rf_ex_phase;
rf_ex.dead_time     = 1e-4;
rf_ex.ringdown_time = 1e-4;
rf_ex.delay         = 1.8e-4;
gz_amplitude        = gz.amplitude;

cd(dir_data)
save([name_rf_exc_pulse],'rf_ex','t_ex','gz_amplitude')

