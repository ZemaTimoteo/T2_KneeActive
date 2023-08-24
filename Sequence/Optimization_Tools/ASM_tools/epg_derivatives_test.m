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


function [ obj, grad, FF ] = epg_derivatives_test( n, tau, R1, R2vec, alpha_RF, alpha )

%% 1 - Initialization
nRates = length(R2vec); % number of T2 to test
tau_half    = tau/2;         % Half of echospacing (ESP/2) in (s)

ds_dT2 = zeros(n,nRates); % number of echoes per number of R2 values
ds_dFl = zeros(n,nRates); % number of echoes per number of R2 values

%% 2 - Selection matrix to move all traverse states up one coherence level
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

%% 3 ==== Cicle for derivatives for points per RF pulse ----
for z=1:size(alpha,1) % |--> pontos dos pulsos (128)
        
    [aux_ds_dT2(:,z),aux_ds_dFl(:,z),aux_epg(:,z)] = epg_derivatives_dEPG( ...
                                    alpha_RF(z), alpha(z), n, nRates, R1, R2vec, tau_half, S ...
                                    );
    
    % 6.3 - ====== Get derivatives for dS/dT2 ==========
%     if nargout > 1 % need gradient
%         
%     end
end

% 3.2 -- Sum of 128 coordinates of z
ds_dT2 = sum(aux_ds_dT2,2);
ds_dFl = sum(aux_ds_dFl,2);
x      = sum(aux_epg,2);

% 3.3 - Plots
% % figure()
% % plot(abs(x))
% % % plot(abs(x./norm(x)))
% % hold on; title('Method toolbox - JUSTdict')

% 3.4 - ====== Outputs ====
obj = 0;                % objective function
grad.ds_dT2 = ds_dT2;   % Gradient for ds/dT2
grad.ds_dFl = ds_dFl;   % Gradient for ds/dFlip_Angle
FF          = x;        % States

