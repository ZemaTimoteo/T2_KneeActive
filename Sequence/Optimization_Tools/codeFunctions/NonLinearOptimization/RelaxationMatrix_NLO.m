function[FpFmZ] = RelaxationMatrix_NLO(FpFmZ,T2,T1,tau)

% 1.2 - Get Relaxation Matrix
% R - relax
    E2 = exp(-tau/T2);
    E1 = exp(-tau/T1);
    EE = diag([E2 E2 E1]);	% Decay of states due to relaxation alone.
    RR = [1-E1];			% Mz Recovery, affects only Z0 state, as
    FpFmZ = EE * FpFmZ;		    % Apply Relaxation
    FpFmZ(3,1) = FpFmZ(3,1)+RR;	% Recovery  ( here applied before diffusion,

end