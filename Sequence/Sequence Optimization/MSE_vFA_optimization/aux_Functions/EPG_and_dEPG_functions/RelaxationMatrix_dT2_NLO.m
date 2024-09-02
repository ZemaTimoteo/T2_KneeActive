function[FpFmZ] = RelaxationMatrix_dT2_NLO(FpFmZ,T2,T1,tau)

% 1.2 - Get Relaxation Matrix
% R - relax
    dT2_E2 = (tau*exp(-tau/T2)) / (T2^2);
    dT2_E1 = 0;
    EE     = diag([dT2_E2 dT2_E2 dT2_E1]);	% Decay of states due to relaxation alone.
    RR     = [1-dT2_E1];                    % Mz Recovery, affects only Z0 state, as
    FpFmZ  = EE * FpFmZ;                    % Apply Relaxation
    FpFmZ(3,1) = FpFmZ(3,1)+RR;             % Recovery  ( here applied before diffusion,

end