function [exc_pulse, refoc_pulse] = slr_profile_fmincon(B1_scale,refoc)
% function [exc_pulse, refoc_pulse] = slr_profile_fmincon(B1_scale,refoc,dTE,acquisition,dir_data)

% Nuno S Santos - 2016
% Excitation and refocusing Shinnar-LeRoux profiles with slice profile correction
% (See paper Buonincontri et al., 2015)


FA_scale = B1_scale;

% Parameters from G and T of excitation pulse and refocusing pulse
G     = 0.74017;  % G/cm - intensity of slice selection -> G = 7.40 mT/m
t     = 2.9376;   % ms - RF pulse duration
G_rf  = 0.61681;  % G/cm - 100 G/cm = 1 T/m = 0.1 G/cm = 1 mT/m -> G_rf = 6.16 mT/m
t_rf  = 1.5232;   % ms
gamma = 42.54;    % Giromagnetic constant (MHz/T)

% Load RF pulses from pypulseq
% load([dir_data '/rf_pulses.mat'])
% puls_exc = rf_ex.signal*gamma;  % RF_ex from pypulseq 

% sg_150_100_167 waveform
% RF_refoc from Phillips
puls_exc = [66, -138, -359, -595, -845, -1107, -1378, -1658, -1943, -2230, -2517, -2800, -3076, -3341, -3592,...
    -3824, -4034, -4217, -4370, -4488, -4567, -4602, -4591, -4530, -4414, -4240, -4007, -3710, -3348, -2918, -2421,...
    -1853, -1216, -510, 266, 1110, 2019, 2991, 4024, 5114, 6256, 7446, 8678, 9947, 11247, 12571, 13911, 15262, 16614,...
    17960,19293, 20604, 21885, 23129, 24326, 25471, 26554, 27570, 28511, 29371,30144, 30826, 31410, 31894, 32274, 32547,...
    32712, 32767, 32712, 32547,32274, 31894, 31410, 30826, 30144, 29371, 28511, 27570, 26554, 25471,24326, 23129, 21885,...
    20604, 19293, 17960, 16614, 15262, 13911, 12571,11247, 9947, 8678, 7446, 6256, 5114, 4024, 2991, 2019, 1110, 266]';

FA_ex    = pi/2*FA_scale;
puls_exc = puls_exc/sum(puls_exc)*FA_ex;

x_ex         = linspace(-5,5,128);
[a_ex, b_ex] = abr(puls_exc,x_ex);

mxy = 2*conj(a_ex).*b_ex;		    % selective excitation
% mz  = 1 - 2*b.*conj(b);	        % inversion
% mxy_se = i*(b.*b);      %spin-echo

flip_angle_z = asin(abs(mxy)); % *(180/pi);


xvec_ex = gt2cm(x_ex,G,t);

%-----------------------------------------------------------------------------------------------------------------
% Echo waveform

% echo 2

% % echo_2 = rf_ref.signal*gamma;      % RF_ex from pypulseq 

% RF_refoc from Phillips
% % echo_2 = [1148, 1605, 1860, 1997, 2070, 2111, 2128, 2118, 2074, 1995, 1879, 1719, 1496, 1186, 764,...
% %     212, -457, -1201, -1960, -2673, -3307, -3828, -4249, -4597, -4902, -5173, -5391, -5507, -5457,...
% %     -5187, -4660, -3869, -2831, -1577, -131, 1491, 3288, 5261, 7465, 9902, 12490, 15180, 17915, ...
% %     20634, 23278, 25778, 28051, 29989, 31490, 32440, 32767, 32440, 31490, 29989, 28051, 25778, 23279,...
% %     20634, 17915, 15181, 12490, 9901, 7465, 5261, 3288, 1491, -130, -1577, -2830, -3868, -4659, -5187,...
% %     -5456, -5506, -5391, -5173, -4902, -4597, -4248, -3827, -3305, -2672, -1960, -1201, -457, 212, 764,...
% %     1187, 1497, 1720, 1880, 1996, 2074, 2118, 2128, 2112, 2071, 1997, 1861, 1606, 1149];

% sinc_centre
echo_2 = [0, 1214, 2536, 3961, 5477, 7074, 8740, 10461, 12223, 14010, 15805, 17592, 19353, 21071, 22728, ...
    24308, 25794, 27170, 28422, 29536, 30500, 31304, 31938, 32397, 32674, 32767, 32674, 32397, 31938, 31304, ...
    30500, 29536, 28422, 27170, 25794, 24308, 22728, 21071, 19353, 17592, 15805, 14010, 12223, 10461, 8740, ...
    7074, 5477, 3961, 2536, 1214, 0];

FA_rf  = refoc*(pi/180)*FA_scale;   % convert refoc angle in rad

echo = echo_2/sum(echo_2)*FA_rf;
x_rf = linspace(-5,5,128);

[~, b_rf] = abr(echo,x_rf);
mxy_se    = 1i*(b_rf.*b_rf);               %spin-echo

flip_angle_z_rf = 2*(asin(sqrt(abs(mxy_se)))); % *(180/pi);

xvec_rf = gt2cm(x_rf,G_rf,t_rf);



% SLICE PROFILE CORRECTION (Buonincontri 2015) ----------------------------------------------------------
% -----------------------------------------------------------------------
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

% ... Interpolate ...
flip_angle_z_cropped_interp = interp1(xvec_longest_cropped,flip_angle_z_cropped,xvec_shortest,'spline'); % n=128 pts


if xvec_ex(end)<xvec_rf(end)
    exc_pulse   = flip_angle_z;
    refoc_pulse = flip_angle_z_cropped_interp;  % (rad)
else
    exc_pulse   = flip_angle_z_cropped_interp;  % (rad)
    refoc_pulse = flip_angle_z_rf;
end

end