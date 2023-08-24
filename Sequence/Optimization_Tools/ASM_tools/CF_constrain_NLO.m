function [x_new,a_k] = CF_constrain_NLO(x_new,params)

% 1 - Parameters - Constrains
constrBetaMin  = params.constr.betaMin;
constrAlphaMin = params.constr.alphaMin;
constrAlphaMax = params.constr.alphaMax;
a_k            = params.constr.t_k;   % itteration factor for guarantee constrain
d_k            = params.constr.d_k;   % gradient value

% 2 - cycle for guarantee constrains
if x_new(1) < constrBetaMin
    x_new(1) = constrBetaMin;
    a_k(1)   = 0;
end
if any(x_new(2:end)< constrAlphaMin)
    idx_cond2        = find(x_new(2:end)< constrAlphaMin == 1);
    x_new(idx_cond2+1) = constrAlphaMin;
    a_k(idx_cond2+1)   = 0;
end
if any(x_new(2:end) > constrAlphaMax)
    idx_cond3        = find(x_new(2:end) > constrAlphaMax == 1);
    x_new(idx_cond3+1) = constrAlphaMax;
    a_k(idx_cond3+1)   = 0;
end


end