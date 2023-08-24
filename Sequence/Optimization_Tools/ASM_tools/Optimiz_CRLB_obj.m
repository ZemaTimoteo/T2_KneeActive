 function [uCRLB] = Optimiz_CRLB_obj(x,constVariables)

%% 0 - parameters 
% variables to optimize
% dTE       = x(1);
% ETL       = x(2);
% flipAngle = x(3);
dTE       = constVariables.dTE;
ETL       = constVariables.ETL;
flipAngle = x(1);

% parameters for EPG
phase_refoc = exp(zeros(1,ETL)./(180/pi).*1i);  % phase of refoc pulses (0 = about x axis)
FA_refoc    = ones(1,ETL)*flipAngle;            % Flip-angle Refocosuing - in (rad) - along y

phase_exc = pi/2;                               % Phase Exciatation - in (rad) - along y
FA_exc    = pi/2;                               % Flip-angle Exciatation - in (rad) - along y

sigma     = (1/constVariables.SNR);              % Noise variance for CRLB


%% 1 - time parameters
% ... 1.1 - Parameters ...
iniTR          = constVariables.TR;    % in (ms)
sliceThickness = constVariables.st;    % in (m)
nsli           = constVariables.nsli;  % number of slices

% ... 1.2 - Get times ...
[TRacq, Trec] = test_CRLB_times_fmincon(dTE,iniTR,ETL,sliceThickness,nsli); % Trec in (ms) & % TRacq in (s)        


%% 2 - Calculate EPG
% ... 2.1 - Parameters for generating signals with slr_profile.py and my_epg.py ...
T1  = constVariables.T1maxKnee;   % ms
T2  = constVariables.T2;          % ms ( 40:1:55; OR 1:1:300; )
B1  = constVariables.B1;          % B1 value ( 0.7:0.01:0.75; OR 0.6:0.01:1.4; )

% ... 2.2 - Obtain dictionary
[S,~] = dict_pars_generator_fmincon(...
                T1,T2,B1,dTE,ETL,...
                phase_refoc,phase_exc,...
                FA_exc,FA_refoc...
                );

% ... 2.3 - Correct Signal ...
aux1 = exp( - Trec/T1);                     % control for Trecovery & TR variation
aux2 = (1 - aux1) / ( 1 - aux1*S(ETL) );
S    = aux2.* S;                            % Signal corrected for variable TR

%% 3 - Calculate uCRLB
    
% ... 3.1 - Parameters ...
step  = 1e-4;                             % step for numerical derivative

% ... 3.2 - Numerical derivatives ...
    % --- 3.2.1 S(M0+step,T2) ---
S_M0_plus_step = (1+step)*S;

    % --- 3.2.2 S(M0,T2+step) ---
refoc_pulse_2 = [];
for jj=1:ETL
    [exc_pulse_2, aux_refoc_pulse_2] = slr_profile_fmincon(B1,flipAngle);
    refoc_pulse_2                    = [refoc_pulse_2 aux_refoc_pulse_2(:)];  
end

S_T2_plus_step  = abs( my_epg_fmincon(exc_pulse_2,refoc_pulse_2,phase_exc,phase_refoc,T1,T2+step,dTE,ETL) );
aux1_S_T2       = exp( - Trec/T1); % control for Trecovery & TR variation
aux2_S_T2       = (1 - aux1_S_T2) / ( 1 - aux1_S_T2 * S_T2_plus_step(ETL) );
S_T2_plus_step  = aux2_S_T2.* S_T2_plus_step;

% ... 3.3 - Num derivatives (forward differentiation) ---
dS_dM0 = (abs(S_M0_plus_step) - abs(S))./(step);
dS_dT2 = (abs(S_T2_plus_step) - abs(S))./(step);

% ... 3.4 - Jacobian Matrix ...
dM_num      = [dS_dM0(:) dS_dT2(:)]; 

% ... 3.5 - CRLB ...
	% --- 3.5.1 CRLBnum ---
FIM_num    = dM_num.'   *  (   (1/sigma^2) * eye( size(dM_num,1) )  )   *  dM_num;
CRLB_num   = diag(  ( inv(FIM_num) )  );

    % --- 3.5.2 Uncertainty of CRLB ---
uCRLB = (abs(CRLB_num(2))/(T2^2))*(TRacq);    % variance CRLB / s - in Zhang et al. 2013 in Proc. EMBS


%% 4 - Get Objective

% % objective = uCRLB;

end