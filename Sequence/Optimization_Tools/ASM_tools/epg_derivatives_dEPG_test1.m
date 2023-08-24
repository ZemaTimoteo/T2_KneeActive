function [ds_dT2,ds_dFl,epg] = epg_derivatives_dEPG_test1(alpha_RF , alpha, n, nRates, R1, R2vec, tau, S)

%% ========================================================================
% Obtain the derivative analytic for epg
%   by: TTFernandes, IST, March 2022
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
%   alpha:  angle of RF (rad)
%   n:      number of echoes
%   nRates: number of T2 tests
%   R1:     1/T1 (s^-1)
%   R2vec:  1/T2 (s^-1)
%   tau:    Echo Time (s)
%   S:      Shift Matrix
%
% Ouputs:
%   ds_dT2: Derivative of ds/dT2
%   ds_dFl: Derivative of ds/dFlip_Angle
%   x:      
%

%% ========================================================================

%% ================== get P(m) ==================================


%% 1 - Get Rotation matrices 
% 1.1 - RF mixing matrix, T (Rotation Matrix)
p_rf = 0;
% % T0 = [ cos(alpha/2)^2, sin(alpha/2)^2,  sin(alpha); ...
% %        sin(alpha/2)^2, cos(alpha/2)^2, -sin(alpha); ...
% %       -0.5*sin(alpha), 0.5*sin(alpha),  cos(alpha)];  % pq é que nao multiplica por 'i', pq é que nao multiplica por exp? Mx, My, Mz
T0 = [ (cos(alpha/2))^2                  ,  exp(2*1i*p_rf)*(sin(alpha/2))^2 , -1i*exp(1i*p_rf)*sin(alpha);  ...
        exp(-2*1i*p_rf)*(sin(alpha/2))^2 , (cos(alpha/2))^2                 ,  1i*exp(-1i*p_rf)*sin(alpha); ...
       -1i/2*exp(-1i*p_rf)*sin(alpha)    ,  1i/2*exp(1i*p_rf)*sin(alpha)    ,  cos(alpha)];

  
% 1.2 - Derivative of T with respect to alpha
% % T0d = [-0.5*sin(alpha),  0.5*sin(alpha),  cos(alpha); ...
% %         0.5*sin(alpha), -0.5*sin(alpha), -cos(alpha); ...
% %        -0.5*cos(alpha),  0.5*cos(alpha), -sin(alpha)];
T0d = [-0.5*sin(alpha)                 ,  0.5*exp(2*1i*p_rf)*sin(alpha) , -1i*exp(1i*p_rf)*cos(alpha); ...
        0.5*exp(-2*1i*p_rf)*sin(alpha) , -0.5*sin(alpha)                ,  1i*exp(-1i*p_rf)*cos(alpha); ...
       -1i/2*exp(-1i*p_rf)*cos(alpha)  ,  1i/2*exp(1i*p_rf)*cos(alpha)  , -sin(alpha)];

% 1.3 - Build mixing matrix for all states
TArray    = cell(1,n);
TArray(:) = {sparse(T0)};
% % T         = blkdiag(1,TArray{:}); % sin alfa
p = pi/2;  
% T         = blkdiag(sin(alpha_RF),TArray{:}); % sin alfa
T         = blkdiag(-1i*exp(1i*p)*sin(alpha_RF),TArray{:}); % sin alfa
% T         = blkdiag(1,TArray{:}); % sin alfa
TArray(:) = {sparse(T0d)};
% Td        = blkdiag(-1i*exp(1i*p)*cos(alpha_RF),TArray{:});
Td        = blkdiag(0,TArray{:});
T_full    = real(full(T));


%% 2 - Loop over different relaxation rates

for iRate=1:nRates % loop over different relaxation rates
    
    % Relaxation matrix
    R2  = R2vec(iRate);
    R0  = diag([exp(-tau*R2),exp(-tau*R2),exp(-tau*R1)]);
    R0d = diag(-tau*exp(-tau*R2)*[1 1 0]);  % derivative w.r.t R2
    
    % Build relaxation matrix (and derivates) for all states
    RArray    = cell(1,n);
    RArray(:) = {sparse(R0)};
    R         = blkdiag(exp(-tau*R2),RArray{:});
    RArray(:) = {sparse(R0d)};
    Rd        = blkdiag(-tau*exp(-tau*R2),RArray{:});
    
    % Precession and relaxation matrix (and derivative)
    P  = (R*S);   % = RS
    Pd = (Rd*S);
    

    % Matrix representing the inter-echo duration
    E    = P*T*P;               % Transition Matrix no derivative
    EdR0 = Pd*T*P + P*T*Pd;     % Transition Matrix with derivatives for T2
    EdR  = EdR0;                % Transition Matrix with derivatives for T2
    EdA0 = P*Td*P;              % Transition matrix for derivatives of flipA
    EdA  = EdA0;                % Transition matrix for derivatives of flipA
    
    % Vector of all states
    x    = zeros(size(R,1),1);
    x(1) = 1;
    
    % full - test 8/4/22
    full_P  = full(P);
    full_T  = full(T);
    full_S  = full(S);
    full_E  = full(E);
    full_R  = full(R);
    full_PT = full(P*T);
    
    full_abs_P  = abs(full_P);
    full_abs_T  = abs(full_T);
    full_abs_S  = abs(full_S);
    full_abs_E  = abs(full_E);
    full_abs_R  = abs(full_R);
    full_abs_PT = abs(full_PT);
    
    % === First echo ===
    % -> Deriviatives of state vector
    xdashR = EdR*x;
    xdashA = EdA*x;
    % -> Update derivative value
    ds_dT2(1,iRate) = xdashR(1);
    ds_dFl(1,iRate) = xdashA(1);
    % -> Calculate new state for next interation
    x = E*x;
    epg(1,iRate) = x(1);
    
    
        % testeee 4/5/22
    x_test    = zeros(size(R,1),1);
    x_test(1) = 1;        
    x_test    = full_T*x_test;
    x_test    = full_R*x_test;
    x_test    = full_S*x_test;
    x_test    = full_T*x_test;
    x_test    = full_R*x_test;
    x_test    = full_S*x_test;
    epg_test(1,iRate) = x_test(1);    
    x = x_test;
    epg(1,iRate) = epg_test(1,iRate);
    
    % === Subsequent echoes ===
    for i=2:n % loop over number of echoes
        % -> Calculate derivatives using the product rule
        xdashR = EdR0*x + E*xdashR; % ( deriv. matrix transit. by states )  +  ( matrix transit. by (deriv. matrix transit. by previous states) )
        xdashA = EdA0*x + E*xdashA;
        % -> Update derivative value
        ds_dT2(i,iRate) = xdashR(1);
        ds_dFl(i,iRate) = xdashA(1);
        % -> Calculate new state for next interation
% %         x = E*x;
        
        % testeee 4/5/22
        x = full_R*x;
        x = full_S*x;
        x = full_T*x;
        x = full_R*x;
        x = full_S*x;
        
        epg(i,iRate) = x(1);
    end
    
end

% % 
% % figure()
% % plot(abs(epg))
end