function [CRLB_num,uCRLB] = CRLB_epg_optFSE(S, phi, refoc_phase, B1, T1, T2, dTE, ETL, flipA, sigma, step, plotTest,dir_data)
% Obtain the variance expected for a specific model/set of parameters
%
% Functions used:
%   slr_profile
%   my_epg: to generate the states over the slice excited (sum) for different T2
%
% Inputs:
%   S:           Original Dictionary signal
%   phi:         phase of excitation
%   refoc_phase: phase of refoc pulses (0 = about x axis)
%   FA_exc:      Flip angle excitation pulse: vector
%   FA_refoc:    Flip angle refocusing pulse: vector
%   B1_scale:    B1 field
%   T1:          constant- relaxation time of the longitudinal magnetization
%   T2:          interval- relaxation time of the transversal magnetization
%   dTE:         Spacing Time between Echoes (ms)
%   ETL:         refocusing echoes train length
%   flipA:       flip angle for refocusing pulse
%   sigma:       Noise variance for CRLB
%   step:        step for numerical derivative
%   plotTest:    'True' or 'Fals'
%
% Ouputs:
%   CRLB_num:    CRLB value
%   uCRLB:       CRLB uncertainty - in Zhang et al. 2013 in Proc. EMBS
%

%%

% ... 1 - Numerical derivatives ...

    % --- 1.1 S(B1+step,T2) ---
refoc_pulse = [];
for jj=1:size(flipA,2)
    [exc_pulse, aux_refoc_pulse] = slr_profile(B1+(step*B1),flipA(jj),dTE,'HSM2',dir_data);
    refoc_pulse                  = [refoc_pulse aux_refoc_pulse(:)];
end

S_B1_plus_step = abs( my_epg(exc_pulse,refoc_pulse,phi,refoc_phase,T1,T2,dTE,ETL) );  


    % --- 1.2 S(B1,T2+step) ---
refoc_pulse_2 = [];
for jj=1:size(flipA,2)
    [exc_pulse_2, aux_refoc_pulse_2] = slr_profile(B1,flipA(jj),dTE,'HSM2',dir_data);
    refoc_pulse_2                    = [refoc_pulse_2 aux_refoc_pulse_2(:)];  
end

S_T2_plus_step = abs( my_epg(exc_pulse_2,refoc_pulse_2,phi,refoc_phase,T1,T2+(step*T2),dTE,ETL) );


    % --- 1.3  Num derivatives (forward differentiation) ---
if plotTest == 'True'
    figure, plot(abs(S),'*--'), hold on, plot(abs(S_B1_plus_step),'ro--'), hold on, plot(abs(S_T2_plus_step),'g+--')
end


    % --- 1.4  Num derivatives (forward differentiation) ---
dS_dB1 = (S_B1_plus_step - S)./(step);
dS_dT2 = (S_T2_plus_step - S)./(step);

% ... 2 - Jacobian Matrix ...
    
dM_num      = [dS_dB1(:) dS_dT2(:)]; 

% ... 3 - CRLB ...

	% --- 3.1 CRLBnum ---
FIM_num    = dM_num.'   *  (   (1/sigma^2) * eye( size(dM_num,1) )  )   *  dM_num;
CRLB_num   = diag(  sqrt( inv(FIM_num) )  );

    % --- 3.2 Uncertainty of CRLB ---
uCRLB = abs(CRLB_num(2))/(T2)/sqrt(nreps);    % in Zhang et al. 2013 in Proc. EMBS


end

