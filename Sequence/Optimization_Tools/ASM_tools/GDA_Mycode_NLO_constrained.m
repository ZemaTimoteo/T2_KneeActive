function [solut] = GDA_Mycode_NLO_constrained(x0,params)

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

% 1.1 - Initialization -----------------------------------------
tol       = params.GDA.tol*ones(size(x0,1));     % tolerence for gradient
a_k       = params.GDA.a_k;     % step of evolution (tip: use a value high enough)
b_k       = params.GDA.b_k;     % step of evolution (tip: use a value high enough)
maxItt    = params.GDA.maxItt;  % Maximum itterations
x_test    = x0;                 % Next point to test - it has dimension alpha_K + 1 (all FA, plus beta)
k         = 1;                  % itterations

% 1.2 - Calculate initial value and initial gradients ----------
[vardT2, aux_g_k] = CF_Mycode_Grad_CRLB_epg_optFSE_NLO(x_test,params);

    % 1.2.1 - Get initial value of cost function
vardT2_test    = vardT2;
vect_vardT2(1) = vardT2_test;
vect_x{1}      = x_test;

    % 1.2.2 - Get initial gradients
g_k          = [aux_g_k.dF_dbeta aux_g_k.dF_dalpha];     % Gradient of [dF_dbeta, dF_dalphaL]
norm_g_k(1)  = norm(g_k);                                % norm of g_k
d_k          = +1 * g_k;                                 % It is a maximization problem, g_k multiplies by (+)

% 1.3 - Summary -------------------------------------------------
fprintf([' -------------------------------\n']);
fprintf(['\n                        Itt: 0 \n\n']);
disp(['Value for xTest of:         - ETL = ', num2str(params.ETL)]);
disp(['                            - TE  = ', num2str(x_test(1)), ' (ms)']);
disp(['Gradient values of  beta is       : ', num2str(g_k(1))]);
disp(['Value for xTest of:         - FA  = ', num2str(x_test(2:end)), ' (radeans)']);
disp(['Gradient values of alpha l is     : ', num2str(g_k(2:end))]);
disp(['Cost Function (vardT2) value is   : ', num2str(vect_vardT2(1))]);
disp(['Gradient (norm g_k) value is      : ', num2str(norm_g_k(1))]);
fprintf([' -------------------------------\n']);
    
    
%% 2 - GDA Loop - Apply Gradient Descend DAlgorithm (GDA)
% 2.0 - Parameters
a_k_test  = a_k;                        % Starting point of parameter evolution (tip: use a value high enough)

while norm_g_k(k) > tol  && k <= maxItt % while gradient isn't smaller than tolerance

    fprintf(['\n\n ---- Test while loop for GDA: ',num2str(k),' itterations ---- \n']);
    
    % 2.1 -  Initial Parameters
    itt        = 0;                      % itterations
    x_new      = x_test + a_k_test.*d_k;  % get x_new
    
    % 2.2 - Garantee constrains ----------------------
    aux_constr = 0;
    aux_itt = 1;
    while aux_constr == 0
        % 2.2.1 - Linear Constrains
        params.constr.t_k = a_k_test;
        params.constr.d_k = d_k;
        [x_new,a_k_test]  = CF_constrain_NLO(x_new,params);    

            % 2.2.2 - Non-linear Constrains
        [params.c,~] = NonLinearConstrains_fmincon_NLOadapted(x_new,params);
        if params.c(1) > 0 && params.c(2) > 0
            aux_constr = 1;
        else
            a_k_test = a_k_test.*b_k;
            x_new    = x_new + a_k_test.*d_k;  % get x_new
        end
        fprintf(['\n Check Constrains Itt:', num2str(aux_itt)]);
        aux_itt = aux_itt+1;
    end
        
    % 2.3 - Get new cost function value --------------
    [vardT2_new,aux_g_k] = CF_Mycode_Grad_CRLB_epg_optFSE_NLO(x_new,params);   % get vardT2_new
        
    % 2.4 - backtracking subroutine ------------------
    while abs(vardT2_new) < vardT2_test && itt < 10 % get funtion of CRLB (x_test + a_k*d_k) > (< - if minimization) function of CRLB (x_test)
        x_new      = x_test + a_k_test.*d_k;
        [vardT2_new,aux_g_k] = CF_Mycode_Grad_CRLB_epg_optFSE_NLO(x_new,params);
        
        % itterate
        a_k_test = a_k_test.*b_k;
        itt = itt + 1; % count itterations
        fprintf(['\n itteration: ',num2str(itt)]);
        
    end
    
    k = k+1;    
   
    
        % ... Check Constrained for initial guess - SAR & Time ...
    [params.c,~] = NonLinearConstrains_fmincon_NLOadapted(x_new,params);
    
    if params.c(1) < 0 || params.c(2) < 0
        fprintf('SAR (uT) check >0:'); display(arams.c(1));
        fprintf('Scan Time (s) check >0:'); display(arams.c(2));
        error('--- Initial guess is out of bounds for SAR or Times (SOLUTION - Design a new initial guess) ---')       
    end
    
    
    
    % 2.5 - Get new values for cost function ---------
    vardT2_test    = vardT2_new;
    vect_vardT2(k) = vardT2_test;
    x_test         = x_new;
    vect_x{k}      = x_test;
    
    % 2.6 - Get gradients -----------------------------
    g_k            = [aux_g_k.dF_dbeta aux_g_k.dF_dalpha];   % Gradient of [dF_dbeta, dF_dalphaL]
    norm_g_k(k)    = norm(g_k);                              % norm of g_k
    d_k            = +1 * g_k;                               % negative of gradient if min, because is max is positive

    % 2.7 - Summary -----------------------------------
    fprintf(['\n                        Itt: ',num2str(itt),'\n\n']);
    disp(['Value for xTest of:         - ETL = ', num2str(params.ETL)]);
    disp(['                            - TE  = ', num2str(x_test(1)), ' (ms)']);
    disp(['Gradient values of  beta is       : ', num2str(g_k(1))]);
    disp(['Value for xTest of:         - FA  = ', num2str(x_test(2:end)), ' (radeans)']);
    disp(['Gradient values of alpha l is     : ', num2str(g_k(2:end))]);
    disp(['Cost Function (vardT2) value is   : ', num2str(vect_vardT2(k))]);
    disp(['Gradient (norm g_k) value is      : ', num2str(norm_g_k(k))]);
    fprintf([' -------------------------------\n']);
        
end

%% 4 - Get Final Solution
solut.x                     = x_test;
solut.vardT2                = vardT2_test;
solut.evolution.vect_vardT2 = vect_vardT2;
solut.evolution.vect_x      = vect_x;
solut.evolution.norm_g_k    = norm_g_k;

end



