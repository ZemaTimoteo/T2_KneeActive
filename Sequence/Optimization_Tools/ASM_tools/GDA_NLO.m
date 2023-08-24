function [sol_vector] = GDA_NLO(x0,params,vardT2)

%% Function to perform Gradient Descent Algorithm (GDA)
% This function will find a minimum - since is a Non Convex problem does
% not garantee a global minimum.
% Will reference it with gradient, and will compare it with the true value
% of the function.
%  Input:
%       - tolerance of gradient     - parms.tol
%       - Starting input for test   - parms.x0
%       - Step                      - parms.step

%% 1 - Implementation

% 1.1 - Initialization
tol       = params.GDA.tol*ones(size(x0,1));     % tolerence for gradient
a_k       = params.GDA.a_k;     % step of evolution (tip: use a value high enough)
b_k       = params.GDA.b_k;     % step of evolution (tip: use a value high enough)
x_test    = x0;                 % Next point to test - it has dimension alpha_K + 1 (all FA, plus beta)
k         = 0;                  % itterations

% 1.2 - Get initial value of cost function
vardT2_test  = vardT2;
   
% 1.3 - Get initial value for gradients
[~, aux_g_k] = CF_Grad_CRLB_epg_optFSE_NLO(x_test,params);
g_k          = [abs(aux_g_k.dF_dbeta) abs(aux_g_k.dF_dalpha)];            % Gradient of [dF_dbeta, dF_dalphaL]
norm_g_k     = norm(g_k);                                       % norm of g_k

d_k = +1 * g_k;                                                 % Como é max, g_k deve ser positivo


fprintf([' -------------------------------\n']);
fprintf(['\n                        Itt: 0 \n\n']);
disp(['Value for xTest of:        - TE = ', num2str(abs(x0(1))), ' (ms)']);
disp(['                           - FA = ' num2str(abs(x0(2:end))), ' (degrees)']);
disp(['Cost Function (vardT2) value is : ' num2str(vardT2_test)]);
disp(['   Gradient (norm g_k) value is : ' num2str(abs(norm_g_k))]);
fprintf([' -------------------------------\n']);
    
    
%% 2 - GDA Loop
% 3 - Apply Gradient Descend DAlgorithm (GDA)
while norm_g_k > tol  % while gradient isn't smaller than tolerance
    
    % 3.1 -  Backtracking routine
    itt        = 0;             % itterations
    a_k_test   = a_k;
    x_new      = x_test + abs(a_k_test*d_k);               % get x_new
    vardT2_new = CF_CRLB_epg_optFSE(x_new,params);    % get vardT2_new  
    
    while abs(vardT2_new) < vardT2_test || itt > 10 % get funtion of CRLB (x_test + a_k*d_k) > (< - if minimization) function of CRLB (x_test)
        x_new      = x_test + a_k_test*d_k;
        vardT2_new = CF_CRLB_epg_optFSE(x_new,params);        
        % itterate
        a_k_test = a_k_test*b_k;       
        itt      = itt + 1; % count itterations
        
    end       
    fprintf(['\n\n Test loop with ',num2str(k),' itterations']);
    
    % 3.2 - set new values for vardT2, x_test, g_k and d_k
    vardT2_test = vardT2_new;     
    x_test      = x_new;
    [~, aux_g_k]    = CF_Grad_CRLB_epg_optFSE_NLO(x_test,params);    % Get Gradient - TODO - colocar funcao de geraçao dos gradientes, nesta linha consoante os parametros de FA e beta(TE)
    g_k         = [aux_g_k.dF_dbeta aux_g_k.dF_dalpha];          % Gradient of [dF_dbeta, dF_dalphaL]
    norm_g_k    = norm(g_k);                                     % norm of g_k
    d_k         = +1 * g_k;                                      % negative of gradient if min, because is max is positive
    
    k = k+1;

    % 3.5 - Summary
fprintf(['\n                        Itt:',num2str(k),'\n\n']);
disp(['Value for xTest of:        - TE = ', num2str(abs(x_test(1))), ' (ms)']);
disp(['                           - FA = ' num2str(abs(x_test(2:end))), ' (degrees)']);
disp(['Cost Function (vardT2) value is : ' num2str(vardT2_test)]);
disp(['   Gradient (norm g_k) value is : ' num2str(abs(norm_g_k))]);
fprintf([' -------------------------------\n']);
    
end



end



