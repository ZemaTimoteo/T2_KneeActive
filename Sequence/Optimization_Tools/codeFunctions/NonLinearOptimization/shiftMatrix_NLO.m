function[FpFmZ] = shiftMatrix_NLO (FpFmZ)

FpFmZ(2,:)   = circshift(FpFmZ(2,:),[0 1]);   % Shift Fm states.
FpFmZ(1,:)   = circshift(FpFmZ(1,:),[0 -1]);  % Shift Fp states.
FpFmZ(1,end) = 0;                             % Zero highest Fp state.
FpFmZ(2,1)   = conj(FpFmZ(1,1));              % Fill in lowest Fm state.
end