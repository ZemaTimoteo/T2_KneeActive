function save_RFref_design(params,st_ref,refoc,dir_data,name_rf_ref_pulse)

%% 1 - Import a seq file to compute SAR for

% ... 1.1 - Parameters for pulseq ...
maxGrad        = 30;                % value for 1.5T & 3T - gm =  30| 32 mT/m,  respectively
maxSlew        = 120;               % value for 1.5T & 3T - sm = 124|130 T/m/s, respectively
sliceThickness = st_ref;            % Slice Thickness (m) - '2.6e-3'
rf_ref_phase   = 0;                 % RF refocusing Phase
fa_val         = refoc;
t_ref          = params.t_ref;      % in ms


% fa_val = [169.4556 162.0682 170.8844 171.3146 168.8563 156.2146 178.991 66.01173 112.3793 143.9194];
% fa_val = [166.8006 139.3508 160.914  179.9874 166.9507 164.4052 147.9974 142.1087 136.22  133.6885 145.7048 157.721 169.7373];
for fa=1:size(fa_val,2)
    % ... 1.2 - initialize sequence ...
    seq = mr.Sequence();
    sys = mr.opts('MaxGrad', maxGrad, 'GradUnit', 'mT/m', 'MaxSlew', maxSlew, ...
        'SlewUnit', 'T/m/s', 'rfRingdownTime', 100e-6, ...
        'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6);
    
    
    %% 2 - Identify RF blocks and compute SAR - 10 seconds must be less than twice and 6 minutes must be less than 4 (WB) and 3.2 (head-20)
    T_vector   = [];
    
    refoc = fa_val(fa);
    % ... 2.2 - Impact of rf_refocusing ...
    % Create rf_refocusing
    flip_ref       = refoc*(pi/180);  % Flip angle refocusing in rad
    [rf_refoc, gz] = mr.makeSincPulse(flip_ref, 'Duration', t_ref,...
        'SliceThickness',sliceThickness,...
        'apodization',0.5, 'timeBwProduct',4,...
        'phaseOffset',rf_ref_phase, 'use','refocusing',...
        'system',sys);
    
    rf_refocDur = mr.calcDuration(rf_refoc);
    signal      = rf_refoc.signal;
    
    %% 3 - Save RF
    % 3.1 - Save RFrefocusing
    rf_ref.type          = 'rf';
    rf_ref.signal        = signal;
    rf_ref.t             = rf_refocDur;
    rf_ref.freq_offset   = 0;
    rf_ref.phase_offset  = rf_ref_phase;
    rf_ref.dead_time     = 1e-4;
    rf_ref.ringdown_time = 1e-4;
    rf_ref.delay         = 2.5e-4;
    rf_ref.use           = 'refocusing';
    gzse_amplitude       = gz.amplitude;
    
    aux_rf_pulse_ref = ['rf_pulses_ref', num2str(refoc)];
    name_rf_pulse    = strrep(aux_rf_pulse_ref,'.','_');

    cd(dir_data)
    save([name_rf_ref_pulse],'rf_ref','t_ref','gzse_amplitude')

end
