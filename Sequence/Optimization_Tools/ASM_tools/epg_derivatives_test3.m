% EPG_DERIVATIVES Derivatives of the Extended Phase Graph (EPG) signal.
%   [Jr Ja] = EPG_DERIVATIVES(N,TAU,R1,R2,ALPHA) calculate the derivatives
%   of the EPG signal with respect to rate (Jr) and flip angle (Ja). The
%   signal is from a CPMG sequence with N echoes and echo spacing of 
%   (2*TAU). The parameters are R1=1/T1 (scalar) and R2=1./T2 (scalar or 
%   vector) and ALPHA the flip angle. The matrices Jr and Ja contains one 
%   column for each element in R2.
%
% The following implementation uses sparse matrices for efficient
% compuation.
%
% Author: Kelvin Layton (klayton@unimelb.edu.au)
% Date: June 2012
% Adapted: TTfernandes (tiagotimoteo@ist.ulisboa.pt)
% Date: March 2022
% Input:
%   - n     : number of Echoes
%   - tau   : echo spacing
%   - R1    : R1 = 1./T1
%   - R2vec : R2 = 1./T2
%   - alpha : flip angle of derivative
%
% Output:
%   - Jr    : Deriviatives of state vector  - ds_dT2
%   - Ja    : Deriviatives of state vector  - ds_dFlip


function [ obj, grad, FF ] = epg_derivatives_test3( n, tau, R1, R2vec, alpha_RF, alpha )

%% 1 - Initialization
nRates = length(R2vec); % number of T2 to test
tau    = tau/2;         % Half of echospacing (ESP/2) in (s)

ds_dT2 = zeros(n,nRates); % number of echoes per number of R2 values
ds_dFl = zeros(n,nRates); % number of echoes per number of R2 values



%% ORIGINAL - Selection matrix to move all traverse states up one coherence level
% Get Shift Matrix
S = sparse(3*n+1,3*n+1);

S(1,3) = 1;
S(2,1) = 1;
S(3,6) = 1;
S(4,4) = 1;

for o=2:n  % loop over number of echoes
    offset1 = ( (o-1) - 1 ) * 3 + 2;
    offset2 = ( (o+1) - 1 ) * 3 + 3;
    if offset1 <= (3*n+1)
        S(3*o-1,offset1) = 1;  % F_k <- F_{k-1}
    end
    if offset2 <= (3*n+1)
        S(3*o,offset2)   = 1;  % F_-k <- F_{-k-1}
    end
    S(3*o+1,3*o+1) = 1;  % Z_order
end

S_fullOriginal = full(S);


%% 2 - Selection matrix to move all traverse states up one coherence level
% Get Shift Matrix
r = 3;
S = sparse(3*n+r,3*n+r);

% S(1,3) = 1;
% S(2,1) = 1;

S(2,5) = 1;
% S(3,5) = 1;
S(4,1) = 1;
S(5,8) = 1;
S(6,6) = 1;

for o=2:n  % loop over number of echoes
    offset1 = ( (o-1) - 1 ) * 3 + 4;
    offset2 = ( (o+1) - 1 ) * 3 + 5;
    if offset1 <= (3*n+3)
        S(3*o+1,offset1) = 1;  % F_k <- F_{k-1}
    end
    if offset2 <= (3*n+3)
        S(3*o+2,offset2)   = 1;  % F_-k <- F_{-k-1}
    end
    S(3*o+3,3*o+3) = 1;  % Z_order
end

S_full3 = full(S);


%% 3 ==== Cicle for derivatives for points per RF pulse ----
for z=1:size(alpha,1) % |--> pontos dos pulsos (128)
        
% %     [aux_ds_dT2(:,z),aux_ds_dFl(:,z),aux_epg] = epg_derivatives_dEPG_test1(alpha_RF(z), alpha(z), n, nRates, R1, R2vec, tau, S);
    [aux_ds_dT2(:,z),aux_ds_dFl(:,z),aux_epg] = epg_derivatives_dEPG_test4(alpha_RF(z), alpha(z), n, nRates, R1, R2vec, tau, S);
    
% %     epg(:,z) = aux_epg./norm(aux_epg);
    s_epg(:,z) = aux_epg;
    % 6.3 - ====== Get derivatives for dS/dT2 ==========
%     if nargout > 1 % need gradient
%         
%     end
end

% 3.2 -- Sum of 128 coordinates of z
ds_dT2 = sum(aux_ds_dT2,2);
ds_dFl = sum(aux_ds_dFl,2);
x      = sum(s_epg,2);

% 6.4 - ====== Outputs ====
obj = 0;                % objective function
grad.ds_dT2 = ds_dT2;   % Gradient for ds/dT2
grad.ds_dFl = ds_dFl;   % Gradient for ds/dFlip_Angle
FF = x;                 % States


% test 29/4
s = s_epg(:,1);
figure; subplot(121); plot(abs(s))
title('Non-norm, Method toolb-SLR, exc=0.0172')
subplot(122); plot(abs(s./norm(s)))
title('Norm, Method toolb.-SLR, exc=0.0172')

aux_s_epg = abs(s_epg);
figure; plot(abs(FF))
