function [exc_pulse, refoc_pulse] = slr_profile_wRFoptimz_design(rf_exc,B1_scale,refoc,params)

% TTFernandes - 2022
% Excitation and refocusing Shinnar-LeRoux profiles with slice profile correction
% (See paper Buonincontri et al., 2015)
% it designs the desired RFexcitation and RF refocus with pulseq design

% Needs pulseq-master path
FA_scale = B1_scale;
dir_data = params.dir_rf;

% Parameters from G and T of excitation pulse and refocusing pulse
t          = params.t_ex;       % ms - RF pulse duration
t_rf       = params.t_ref;      % ms
gamma      = 42.54;             % Giromagnetic constant (MHz/T)
center_pos = params.center_pos; % center_pos
npoints    = params.npoints;

%% ------------------------------------------------------------------------
% 1 - Load OR Create RF excitation and refocusing
cd(dir_data);

% 1.1 - Load OR Create RFexcit
if params.RFqual == 'True'
    aux_name_rf_exc_pulse = ['rf_pulses_exc', num2str(rf_exc*180/pi), '_assym',num2str((center_pos-0.5)*100),'perCent_st', num2str(params.st_ex*1e3),'mm'];
    st_ex = params.st_ex;   
else
    aux_name_rf_exc_pulse = ['rf_pulses_exc', num2str(rf_exc*90/(pi/2)), '_assym',num2str((center_pos-0.5)*100),'perCent'];
    st_ex = params.st;    
end
aux_name_rf_exc_pulse = strrep(aux_name_rf_exc_pulse,'.','_');
name_rf_exc_pulse     = [aux_name_rf_exc_pulse '.mat'];
if isfile(name_rf_exc_pulse)
    load(name_rf_exc_pulse);
else
    save_RFexc_design(params,st_ex,rf_exc,dir_data,name_rf_exc_pulse);
    load(name_rf_exc_pulse);
end

% 1.2 - Load OR Create RFrefoc
if params.RFqual == 'True'
    aux_name_rf_ref_pulse = ['rf_pulses_ref', num2str(refoc), '_st', num2str(params.st_ref*1e3),'mm']; 
    st_ref = params.st_ref;
else
    aux_name_rf_ref_pulse = ['rf_pulses_ref', num2str(refoc)];
    st_ref = params.st;    
end
aux_name_rf_ref_pulse = strrep(aux_name_rf_ref_pulse,'.','_');
name_rf_ref_pulse = [aux_name_rf_ref_pulse '.mat'];
if isfile(name_rf_ref_pulse)
    load(name_rf_ref_pulse);
else
    save_RFref_design(params,st_ref,refoc,dir_data,name_rf_ref_pulse);
    load(name_rf_ref_pulse);
end


%% ------------------------------------------------------------------------
% 2 - RF excitation waveform
    % 2.1 - Load RF excitation pulses from pypulseq
puls_exc = rf_ex.signal*gamma;  % RF_ex from pypulseq 
FA_ex    = rf_exc*FA_scale;
puls_exc = puls_exc/sum(puls_exc)*FA_ex;
G_ex     = gz_amplitude/(gamma*1e3) * 1e-1;           % G/cm - intensity of slice selection -> 0.74017 -> G = 7.40 mT/m
% puls_exc = rf_ex.signal*gamma * 1e-3 * FA_scale;      % mT

x_ex         = linspace(-5,5,npoints);
[a_ex, b_ex] = abr(puls_exc,x_ex);

mxy          = 2*conj(a_ex).*b_ex;	 % selective excitation
flip_angle_z = asin(abs(mxy));       % *(180/pi);

xvec_ex = gt2cm(x_ex,G_ex,t);

%% ------------------------------------------------------------------------
% 3 - Echo waveform

    % 3.1 - Load RF refoc pulse from pypulseq
echo_2 = rf_ref.signal*gamma;  % RF_refoc from pypulseq 
FA_rf  = refoc*(pi/180)*FA_scale;   % convert refoc angle in rad
echo   = echo_2/sum(echo_2)*FA_rf;
G_rf   = gzse_amplitude/(gamma*1e3) * 1e-1;         % G/cm - intensity of slice selection -> 0.74017 -> G = 7.40 mT/m
% echo   = rf_ref.signal*gamma * 1e-3 *  FA_scale;     % mT

x_rf      = linspace(-5,5,npoints);
[~, b_rf] = abr(echo,x_rf);

mxy_se          = 1i*(b_rf.*b_rf);               %spin-echo
flip_angle_z_rf = 2*(asin(sqrt(abs(mxy_se))));   % *(180/pi);

xvec_rf = gt2cm(x_rf,G_rf,t_rf);

%% ------------------------------------------------------------------------
% 4 - SLICE PROFILE CORRECTION (Buonincontri 2015) ------------------------
% check the longest profile to crop and interpolate (Exc or Refoc)

if xvec_ex(end)<xvec_rf(end)
    xvec_longest       = xvec_rf;
    xvec_shortest      = xvec_ex;
    flip_angle_longest = flip_angle_z_rf;
else
    xvec_longest       = xvec_ex;
    xvec_shortest      = xvec_rf;
    flip_angle_longest = flip_angle_z;
end

min_diff_profiles = min(abs(abs(xvec_longest) - xvec_shortest(end)));
if min_diff_profiles>1e-2
    error(['pulses are very different - ' num2str(min_diff_profiles) ', and Npoints=' num2str(npoints)])
end

ind  = find(abs(abs(xvec_longest) - xvec_shortest(end))<1e-2);
dif  = abs(abs(xvec_longest(ind)) - xvec_shortest(end));
ind2 = find(dif==max(dif));  % min em vez de maximo
    
xvec_longest_cropped = xvec_longest(ind(ind2(1)):ind(ind2(end)));    
flip_angle_z_cropped = flip_angle_longest(ind(ind2(1)):ind(ind2(end)));  % n=98 pts

%% ------------------------------------------------------------------------
% 4 - Interpolate
flip_angle_z_cropped_interp = interp1(xvec_longest_cropped,flip_angle_z_cropped,xvec_shortest,'spline'); % n=128 pts

if xvec_ex(end)<xvec_rf(end)
    exc_pulse   = flip_angle_z;
    refoc_pulse = flip_angle_z_cropped_interp;  % (rad)
else
    exc_pulse   = flip_angle_z_cropped_interp;  % (rad)
    refoc_pulse = flip_angle_z_rf;
end

end