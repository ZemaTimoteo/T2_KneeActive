function [F] = auxOptimiz_epg_cpmg(np, alph, ph, RS, N)

%% ========================================================================
% Obtain the epg_cpmg for excitation and refocusing pulses in MSE
%   by: TTFernandes, IST, Feb. 2022
%
% Functions used:
%   Trot_fun_ASM.m
%   blockdiag_mult_ASM.m
%
% Inputs:
%   np:   number of states
%   alph: Angle vector
%   ph:   Phase vector
%   RS:   composites of R & S matrices
%   N:    number of states in total - F+; F-; Mz p/echo
%
% Ouputs:
%   F:    CRLB value
%

%% ========================================================================

%% 1 - Initialization
F      = zeros([N*2 np+1]);  % records state
F(3,1) = 1;                  % initial state:

for jj = 2:np+1 %loop over time
    
    % excitation
    if jj == 2
        % Excitation pulse
        A                   = Trot_fun_ASM(alph(jj),ph(jj)); %<---- Rotation matrix direct definition
        F([1:3 N+1:N+3],jj) = [real(A); imag(A)] * F(1:3,jj-1);    %<---- state straight after excitation [F0 F0* Z0]
        
        % refocusing
    else
        % Relaxation before RFocus
        temp1    = RS * F(1:N,jj-1);
        temp2    = RS * F(N+1:2*N,jj-1);
        temp2(1) = -temp2(1);% D*temp;
        
        % Rotation with RFocus
        A = Trot_fun_ASM(alph(jj),ph(jj)); % Rotation matrix direct definition
        F(:,jj)  = [ blockdiag_mult_ASM(real(A),temp1) - blockdiag_mult_ASM(imag(A),temp2); ...
            blockdiag_mult_ASM(imag(A),temp1) + blockdiag_mult_ASM(real(A),temp2) ];
        
        % Relaxation after RFocus
        F(1:N,jj)     = RS * F(1:N,jj);
        F(N+1:2*N,jj) = RS * F(N+1:2*N,jj);

    end
end

end