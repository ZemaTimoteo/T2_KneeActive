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
%   - n        : number of Echoes
%   - tau      : echo spacing
%   - R1       : R1 = 1./T1
%   - R2vec    : R2 = 1./T2
%   - alpha_RF : excitation RF pulse
%   - alpha    : refocusing flip angle of derivative
%
% Output:
%   - Jr    : Deriviatives of state vector  - ds_dT2
%   - Ja    : Deriviatives of state vector  - ds_dFlip


function [ grad ] = CF_epg_derivatives_NLO( ETL, tau, R1, R2vec, alpha_RF, alpha, theta)
%% 1 - Initialization
nRates = length(R2vec); % number of T2 to test
tau    = tau/2;         % Half of echospacing (ESP/2) in (s) = TE/2
ds_dT2 = zeros(ETL,nRates); % number of echoes per number of R2 values

%% 2 - Selection matrix to move all traverse states up one coherence level
% Get Shift Matrix
S = sparse(3*ETL+1,3*ETL+1);

S(1,3) = 1;
S(2,1) = 1;
S(3,6) = 1;
S(4,4) = 1;

for o=2:ETL  % loop over number of echoes
    offset1 = ( (o-1) - 1 ) * 3 + 2;
    offset2 = ( (o+1) - 1 ) * 3 + 3;
    if offset1 <= (3*ETL+1)
        S(3*o-1,offset1) = 1;  % F_k <- F_{k-1}
    end
    if offset2 <= (3*ETL+1)
        S(3*o,offset2)   = 1;  % F_-k <- F_{-k-1}
    end
    S(3*o+1,3*o+1) = 1;  % Z_order
end

S_full = full(S);

%% 3 ==== Cicle for derivatives for points per RF pulse ----
for z=1:size(alpha,1) % |--> pontos dos pulsos (128)
        
    [aux_ds_dT2(:,z),aux_dF_dalpha(:,:,z),aux_dF_dbeta(:,z)] = CF_epg_derivatives_aux_NLO( ...
                                        alpha_RF(z,:), alpha(z,:), ETL, nRates, R1, ...
                                        R2vec, tau, S, theta);
end

% 3.2 -- Sum of 128 coordinates of z
ds_dT2    = sum(aux_ds_dT2,2);
dF_dalpha = sum(aux_dF_dalpha,3);
dF_dbeta  = sum(aux_dF_dbeta,2);

% 6.4 - ====== Outputs ====
grad.ds_dT2    = ds_dT2;    % Gradient for ds/dT2
grad.dF_dalpha = dF_dalpha; % Gradient for dF_dalpha
grad.dF_dbeta  = dF_dbeta;  % Gradient for dF_dbeta

