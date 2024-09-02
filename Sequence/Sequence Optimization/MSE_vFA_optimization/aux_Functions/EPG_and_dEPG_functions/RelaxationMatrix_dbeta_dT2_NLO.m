function[FpFmZ] = RelaxationMatrix_dbeta_dT2_NLO(FpFmZ,T2,T1,tau)

% 1.2 - Get Relaxation Matrix
% R - relax
    dbeta_dT2_E2 =  ( exp(-tau/T2) * (T2 - tau) )  / (T2^3);
    dbeta_dT2_E1 = 0;
    EE           = diag([dbeta_dT2_E2 dbeta_dT2_E2 dbeta_dT2_E1]);	% Decay of states due to relaxation alone.
    RR           = [1-dbeta_dT2_E1];                                % Mz Recovery, affects only Z0 state, as
    FpFmZ        = EE * FpFmZ;                                      % Apply Relaxation
    FpFmZ(3,1)   = FpFmZ(3,1)+RR;                                   % Recovery  ( here applied before diffusion,

end