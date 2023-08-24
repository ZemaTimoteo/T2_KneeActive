function[FpFmZ] = RelaxationMatrix_dBeta_NLO(FpFmZ,T2,T1,tau)

% 1.2 - Get Relaxation Matrix
% R - relax
    dbeta_E2   = -exp(-tau/T2) / T2;
    dbeta_E1   = -exp(-tau/T1) / T1;
    EE         = diag([dbeta_E2 dbeta_E2 dbeta_E1]);	% Decay of states due to relaxation alone.
    RR         = [1-dbeta_E1];                          % Mz Recovery, affects only Z0 state, as
    FpFmZ      = EE * FpFmZ;                            % Apply Relaxation
    FpFmZ(3,1) = FpFmZ(3,1)+RR;                         % Recovery  ( here applied before diffusion,

end