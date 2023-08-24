




function [Matrix_dS_dT2] = auxOptimiz_dEPG(alph, ph, RS_full, T2, TE, S0, N, S, np, FF)

%% ========================================================================
% Obtain the derivative analytic for epg
%   by: TTFernandes, IST, Feb. 2022
%
%%% ================== Math ==============================================
% S  = SUM( Fn * Cn )
% Fn = P(n-1) * F(n-1)
% P  = I (x) (RE)S
%
% dS/dT2 = SUM{n=1 to N-1}     [ PI{m=N-1 to n+1}[ P(m) ] ]    *    [ dP/dT2(n) ]   *   [F(n)]
%%% ======================================================================
%
%
% Functions used:
%   Trot_fun_ASM.m
%   blockdiag_mult_ASM.m
%   dErelax_fun_dT2_ASM.m
%
% Inputs:
%   np:   number of states
%   ns:   loop over space
%   alph: Angle vector
%   ph:   Phase vector
%   RS:   composites of R & S matrices
%   N:    number of states in total - F+; F-; Mz p/echo
%
% Ouputs:
%   F:    CRLB value
%

%% ========================================================================

%% ================== get P(m) ===================================

for jj = 1:np+1 %loop over time    
    A      = Trot_fun_ASM(alph(jj),ph(jj)); % Rotation matrix direct definition
    if jj ==1      % Zero state
        P{jj} = zeros(3,3);
    elseif jj == 2 % Excitation %<--- P  = R
        P{jj} = A;
    else           % Refocusing %<--- P  = I (x) (RE)S
        P{jj}  = A*RS_full(1:3,1:3);
    end
end

%% ================== get dP/dT2(n) =============================
% get derivative relaxation matrix w.r.t T2
dE_dT2eff      = dErelax_fun_dT2_ASM(T2,TE/2,S0,N); % derivative w.r.t T2
dRS            = dE_dT2eff*S;

dE_dT2eff_full = dErelax_fun_dT2_ASM(T2,TE,S0,N); % derivative w.r.t T2
dRS_full       = dE_dT2eff_full*S;

for jj = 1:np+1 %loop over time    
    A          = Trot_fun_ASM(alph(jj),ph(jj)); % Rotation matrix direct definition
    if jj ==1      % Zero state
        dP_dT2{jj} = zeros(3,3);
    elseif jj == 2 % derivative of Excitation     %<--- dP/dT2  = dR/dT2 = 0
        dP_dT2{jj} = zeros(3,3);
% %     elseif jj == 3 % derivative of 1st Refocusing %<--- dP/dT2  = I (x) (R * dE/dT2)S for TE/2
% %         dP_dT2{jj} = A*dRS(1:3,1:3);
    else           % derivative of Refocusing     %<--- dP/dT2  = I (x) (R * dE/dT2)S for TE
        dP_dT2{jj} = A*dRS_full(1:3,1:3);
    end
end

clear dE_dT2eff dE_dT2eff_full


%% ================== get dS/dT2(n) with all states ===========
%  get cycle of m variable 
%  get cycle of n variable 
% dS/dT2 = SUM{n=1 to N-1}     [ PI{m=N-1 to n+1}[ P(m) ] ]    *    [ dP/dT2(n) ]   *   [F(n)]


% initialize
nn = N;
Matrix_dS_dT2 = zeros(nn,np+1);

for N = 2:np+2 % loop over number of different states
    idx_N = N-1;
    
    if N==2
        Matrix_dS_dT2(:,N) = blockdiag_mult_ASM(dP_dT2{N-1},FF(1:nn,N-1)) ;
        
    else
        aux_dS_dT2    = zeros(nn,1);
        for n = N-2:-1: 1 % loop within derivative
            % --- Set 1st term of expression ---
            clear P_aux
            P_aux = 1;
            for m=N-1:-1:n+1
                P_aux = P_aux*P{m};
            end
            
            % --- Set 2nd term of expression ---
            aux_Matrix = blockdiag_mult_ASM(P_aux,blockdiag_mult_ASM(dP_dT2{n},FF(1:nn,n)));
            aux_dS_dT2 = aux_Matrix + aux_dS_dT2 ;
            
            clear P_aux
        end
        
        % Add P(m) independent part
        Matrix_dS_dT2(:,N) = blockdiag_mult_ASM(dP_dT2{N-1},FF(1:nn,N-1)) + aux_dS_dT2 ;
        test_dS_dT2(:,N)   = aux_dS_dT2;
        clear aux_dS_dT2
    end
end

end