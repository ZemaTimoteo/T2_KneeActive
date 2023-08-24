function [exc_pulse, refoc_pulse] = slr_profile_vF(rf_exc,B1_scale,refoc,dir_data)

% TTFernandes - 2022
% Excitation and refocusing Shinnar-LeRoux profiles with slice profile correction
% (See paper Buonincontri et al., 2015)


FA_scale = B1_scale;

% Parameters from G and T of excitation pulse and refocusing pulse
G     = 0.74017;  % G/cm - intensity of slice selection -> G = 7.40 mT/m
t     = 2.9376;   % ms - RF pulse duration
G_rf  = 0.61681;  % G/cm - 100 G/cm = 1 T/m = 0.1 G/cm = 1 mT/m -> G_rf = 6.16 mT/m
t_rf  = 1.5232;   % ms
gamma = 42.54;    % Giromagnetic constant (MHz/T)

%% ------------------------------------------------------------------------
% 1 - RF excitation waveform
    % 1.1 - Load RF excitation pulses from pypulseq
dir_rf_exc_pulse = [dir_data,'\rf_pulses_exc', num2str(rf_exc*180/pi),'.mat']; load(dir_rf_exc_pulse);
puls_exc         = rf_ex.signal*gamma;  % RF_ex from pypulseq 

FA_ex    = pi/2*FA_scale;
puls_exc = puls_exc/sum(puls_exc)*FA_ex;

x_ex         = linspace(-5,5,128);
[a_ex, b_ex] = abr(puls_exc,x_ex);

mxy          = 2*conj(a_ex).*b_ex;	 % selective excitation
flip_angle_z = asin(abs(mxy));       % *(180/pi);

xvec_ex = gt2cm(x_ex,G,t);

%% ------------------------------------------------------------------------
% 2 - Echo waveform

    % 2.1 - Load RF refoc pulse from pypulseq
dir_rf_refoc_pulse = [dir_data,'\rf_pulses_ref', num2str(refoc),'.mat']; load(dir_rf_refoc_pulse);
echo_2             = rf_ref.signal*gamma;  % RF_refoc from pypulseq 

FA_rf  = refoc*(pi/180)*FA_scale;   % convert refoc angle in rad
echo   = echo_2/sum(echo_2)*FA_rf;

x_rf      = linspace(-5,5,128);
[~, b_rf] = abr(echo,x_rf);

mxy_se          = 1i*(b_rf.*b_rf);               %spin-echo
flip_angle_z_rf = 2*(asin(sqrt(abs(mxy_se))));   % *(180/pi);

xvec_rf = gt2cm(x_rf,G_rf,t_rf);

%% ------------------------------------------------------------------------
% 3 - SLICE PROFILE CORRECTION (Buonincontri 2015) ------------------------
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