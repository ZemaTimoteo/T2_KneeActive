 function [ds_dT2,ds_dFl,epg,grad] = epg_derivatives_fmincon_dEPG_vF(alpha_RF , alpha, n, nRates, R1, R2vec, tau, S)             

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
%   alpha:  angle of RF (FA) (rad)
%   n:      number of echoes
%   nRates: number of T2 tests
%   R1:     1/T1 (s^-1)
%   R2vec:  1/T2 (s^-1)
%   tau:    Echo Time (s) - TE
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
p_ex = pi/2;  
T0   =  [ (cos(alpha/2))^2                  ,  exp(2*1i*p_rf)*(sin(alpha/2))^2 , -1i*exp(1i*p_rf)*sin(alpha);  ...
           exp(-2*1i*p_rf)*(sin(alpha/2))^2 , (cos(alpha/2))^2                 ,  1i*exp(-1i*p_rf)*sin(alpha); ...
          -1i/2*exp(-1i*p_rf)*sin(alpha)    ,  1i/2*exp(1i*p_rf)*sin(alpha)    ,  cos(alpha)];
TArray    = cell(1,n);
TArray(:) = {sparse(T0)};
T         = blkdiag(-1i*exp(1i*p_ex)*sin(alpha_RF),TArray{:}); % sin alfa
  
% 1.2 - Derivative of T with respect to alpha
T0d = [-0.5*sin(alpha)                 ,  0.5*exp(2*1i*p_rf)*sin(alpha) , -1i*exp(1i*p_rf)*cos(alpha); ...
        0.5*exp(-2*1i*p_rf)*sin(alpha) , -0.5*sin(alpha)                ,  1i*exp(-1i*p_rf)*cos(alpha); ...
       -1i/2*exp(-1i*p_rf)*cos(alpha)  ,  1i/2*exp(1i*p_rf)*cos(alpha)  , -sin(alpha)];

% 1.3 - Build mixing matrix for all states
TArray(:) = {sparse(T0d)};
Td        = blkdiag(0,TArray{:});
T_full    = real(full(T));

       % --- Matrix derivatives for fmincon ---
% 1.1 - RF mixing matrix, T (Rotation Matrix) derivate to FA
T0_dFA        =  T0d;
TArray_dFA    = cell(1,n);
TArray_dFA(:) = {sparse(T0_dFA)};
T_dFA         = blkdiag(0,TArray_dFA{:}); % sin alfa

% 1.1 - RF mixing matrix, T (Rotation Matrix) derivate to TE
T0_dTE        = zeros(3,3);
TArray_dTE    = cell(1,n);
TArray_dTE(:) = {sparse(T0_dTE)};
T_dTE         = blkdiag(0,TArray_dTE{:});

% 1.1 - RF mixing matrix, T (Rotation Matrix) derivate to ETL
T0_dETL        = zeros(3,3);
TArray_dETL    = cell(1,n);
TArray_dETL(:) = {sparse(T0_dETL)};
T_dETL         = blkdiag(0,TArray_dETL{:});

% 1.2 - Derivative of T with respect to alpha
T0d_dFA = [-0.5*cos(alpha)                 ,  0.5*exp(2*1i*p_rf)*cos(alpha) ,   1i*exp(1i*p_rf)*sin(alpha); ...
            0.5*exp(-2*1i*p_rf)*cos(alpha) , -0.5*cos(alpha)                ,  -1i*exp(-1i*p_rf)*sin(alpha); ...
            1i/2*exp(-1i*p_rf)*sin(alpha)  , -1i/2*exp(1i*p_rf)*sin(alpha)  , -cos(alpha)];
T0d_dTE  = zeros(3,3);
T0d_dETL = zeros(3,3);
   
% 1.3 - Build mixing matrix for all states
TArray_dFA(:) = {sparse(T0d_dFA)};
Td_dFA        = blkdiag(0,TArray_dFA{:});
T_full_dFA    = real(full(T_dFA));

TArray_dTE(:) = {sparse(T0d_dTE)};
Td_dTE        = blkdiag(0,TArray_dTE{:});
T_full_dTE    = real(full(T_dTE));

TArray_dETL(:) = {sparse(T0d_dETL)};
Td_dETL        = blkdiag(0,TArray_dETL{:});
T_full_dETL    = real(full(T_dETL));
       % ---------------------------------------

       
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
    
    
           % --- Matrix derivatives for fmincon ---
    % Relaxation matrix
    R2_dFA  = R2vec(iRate);
    R0_dFA  = diag([0,0,0]);
    R0d_dFA = R0_dFA;  % derivative w.r.t R2

    R2_dTE  = R2vec(iRate);
    R0_dTE  = diag([-R2*exp(-tau*R2),-R2*exp(-tau*R2),-R1*exp(-tau*R1)]);
    R0d_dTE = R0_dTE;  % derivative w.r.t R2
    
    R2_dETL  = R2vec(iRate);
    R0_dETL  = diag([0,0,0]);
    R0d_dETL = R0_dETL;  % derivative w.r.t R2
    
    % Build relaxation matrix (and derivates) for all states
    RArray_dFA    = cell(1,n);
    RArray_dFA(:) = {sparse(R0_dFA)};
    R_dFA         = blkdiag(0,RArray_dFA{:});
    RArray_dFA(:) = {sparse(R0d_dFA)};
    Rd_dFA        = blkdiag(0,RArray_dFA{:});

    RArray_dTE    = cell(1,n);
    RArray_dTE(:) = {sparse(R0_dTE)};
    R_dTE         = blkdiag(0,RArray_dTE{:});
    RArray_dTE(:) = {sparse(R0d_dTE)};
    Rd_dTE        = blkdiag(0,RArray_dTE{:});
    
    RArray_dETL    = cell(1,n);
    RArray_dETL(:) = {sparse(R0_dETL)};
    R_dETL         = blkdiag(0,RArray_dETL{:});
    RArray_dETL(:) = {sparse(R0d_dETL)};
    Rd_dETL        = blkdiag(0,RArray_dETL{:});    
    
    % Precession and relaxation matrix (and derivative)
    %     P  = (R*S);   % = RS
    P_dFA  = R_dFA*S;
    P_dTE  = R_dTE*S;
    P_dETL = R_dETL*S;
    
    %     Pd = (Rd*S);
    Pd_dFA  = Rd_dFA*S;
    Pd_dTE  = Rd_dTE*S;
    Pd_dETL = Rd_dFA*S;

           % ---------------------------------------

           
    % Matrix representing the inter-echo duration
    E    = P*T*P;               % Transition Matrix no derivative
    
    % Transition Matrix with derivatives for T2
    EdR0 = Pd*T*P + P*T*Pd;     % Transition Matrix with derivatives for T2
    EdR  = EdR0;                % Transition Matrix with derivatives for T2    
    
    % Transition matrix for derivatives of flipA
    EdA0 = P*Td*P;              % Transition matrix for derivatives of flipA
    EdA  = EdA0;                % Transition matrix for derivatives of flipA
    
    
       % --- Matrix derivatives for fmincon ---
    % Matrix representing the inter-echo duration
       %    E    = P*T*P;
    E_dFA  = P_dFA*T*P + P*T_dFA*P + P*T*P_dFA;               % Transition Matrix derivative         
    E_dTE  = P_dTE*T*P + P*T_dTE*P + P*T*P_dTE;               % Transition Matrix derivative         
    E_dETL = P_dFA*T*P + P*T_dFA*P + P*T*P_dFA;               % Transition Matrix derivative         

    % Transition Matrix with derivatives for T2
       %    EdR0 = Pd*T*P + P*T*Pd;
    EdR0_dFA  = Pd_dFA*T*P + Pd*T_dFA*P + Pd*T*P_dFA   + ...
                   P_dFA*T*Pd + P*T_dFA*Pd + P*T*Pd_dFA;
    EdR0_dTE  = Pd_dTE*T*P + Pd*T_dTE*P + Pd*T*P_dTE   + ...
                   P_dTE*T*Pd + P*T_dTE*Pd + P*T*Pd_dTE;
    EdR0_dETL = Pd_dETL*T*P + Pd*T_dETL*P + Pd*T*P_dETL   + ...
                   P_dETL*T*Pd + P*T_dETL*Pd + P*T*Pd_dETL;               
       % ---------------------------------------
       
       
    % Vector of all states
    x    = zeros(size(R,1),1);
    x(1) = 1;
    
    % === RF Excitation ===
    x = T*x;    
    
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

    
        % --- Matrix derivatives for fmincon ---    
    x_dFA  = zeros(size(R,1),1);
    x_dTE  = zeros(size(R,1),1);
    x_dETL = zeros(size(R,1),1);
    
    x_dFA(1)  = 0;  
    x_dTE(1)  = 0;  
    x_dETL(1) = 0;  
    
    % === RF Excitation ===
    %     x = T*x;    
    x_dFA  = T_dFA*x  + T*x_dFA;
    x_dTE  = T_dTE*x  + T*x_dTE;
    x_dETL = T_dETL*x + T*x_dETL;
    
    % === First echo ===   
    % -> Deriviatives of state vector
    %     x = E*x;    
    xdashR_dFA  = EdR0_dFA*x  + EdR0*x_dFA;
    xdashR_dTE  = EdR0_dTE*x  + EdR0*x_dTE;
    xdashR_dETL = EdR0_dETL*x + EdR0*x_dETL;
    
    % -> Calculate new derivate state for next interation
    x_dFA  = E_dFA*x  + E*x_dFA;
    x_dFA  = E_dFA*x  + E*x_dFA;
    x_dETL = E_dETL*x + E*x_dETL;
       % -----------------------------------------------

    
    % === Subsequent echoes ===
    for i=2:n % loop over number of echoes
        
        % -> Calculate derivatives using the product rule
        xdashR = EdR0*x + E*xdashR; % ( deriv. matrix transit. by states )  +  ( matrix transit. by (deriv. matrix transit. by previous states) )
        xdashA = EdA0*x + E*xdashA;
        % -> Update derivative value
        ds_dT2(i,iRate) = xdashR(1);
        ds_dFl(i,iRate) = xdashA(1);
        % -> Calculate new state for next interation
        x = E*x;  
        % -> Get EPG        
        epg(i,iRate) = x(1);
        
        
            % --- Matrix derivatives for fmincon ---    
        % -> Calculate derivatives for x
        %   x = E*x;
        x_dFA  = E_dFA*x  + E*x_dFA;
        x_dTE  = E_dTE*x  + E*x_dTE;
        x_dETL = E_dETL*x + E*x_dETL;
                    
        % -> Calculate derivatives using the product rule
        %    xdashR = EdR0*x + E*xdashR;        
        xdashR_dFA  = EdR0_dFA*x + EdR0*x_dFA    + ...
                      E_dFA*xdashR + E*xdashR_dFA;
        xdashR_dTE  = EdR0_dTE*x + EdR0*x_dTE    + ...
                      E_dTE*xdashR + E*xdashR_dTE;                  
        xdashR_dETL = EdR0_dETL*x + EdR0*x_dETL  + ...
                      E_dETL*xdashR + E*xdashR_dETL;
                  
        % -> Update Gradient for objetive function   
        %    ds_dT2(i,iRate) = xdashR(1);
        grad.ds_dT2_dFA(i,iRate)  = xdashR_dFA(1);
        grad.ds_dT2_dTE(i,iRate)  = xdashR_dTE(1);
        grad.ds_dT2_dETL(i,iRate) = xdashR_dETL(1);

        
            % ---------------------------------------------        
        
    end
    
end

end