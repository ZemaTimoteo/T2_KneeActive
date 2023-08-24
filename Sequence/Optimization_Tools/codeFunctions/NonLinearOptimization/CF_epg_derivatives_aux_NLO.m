 function [ds_dT2,df_dalpha,df_dbeta] = CF_epg_derivatives_aux_NLO(alpha_RF , alpha, n, nRates, R1, R2vec, tau, S, theta)             

%% ========================================================================
% Obtain the derivative analytic for epg
%   by: TTFernandes, IST, June 2023
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

%% 1 - Get Rotation matrices (T) and dT_partial_alpha

p_rf = 0;
p_ex = pi/2;  

% 1.1 - RF mixing matrix, T (Rotation Matrix)
for j=1:n % loop over number of echoes 
    Tzeros0 = zeros(3,3);
    T0   =  [ (cos(alpha(j)/2))^2                  ,  exp(2*1i*p_rf)*(sin(alpha(j)/2))^2 , -1i*exp(1i*p_rf)*sin(alpha(j));  ...
               exp(-2*1i*p_rf)*(sin(alpha(j)/2))^2 , (cos(alpha(j)/2))^2                 ,  1i*exp(-1i*p_rf)*sin(alpha(j)); ...
              -1i/2*exp(-1i*p_rf)*sin(alpha(j))    ,  1i/2*exp(1i*p_rf)*sin(alpha(j))    ,  cos(alpha(j))];
    % Ge in cell arrays   
    TArray{j}      = T0;
    TzerosArray{j} = Tzeros0;
end

% 1.2 - Derivative of T with respect to alpha
namesGenerate         = 1:n;
TdArray_partial_alpha = [];
Td_partial_alpha      = [];

for j=1:n % loop over number of echoes               
    % Initialization
%     T0partialFAdArray = setfield(T0partialFAdArray, ['fa_',num2str(namesGenerate(j))], TArray);
    TdArray_partial_alpha = setfield(TdArray_partial_alpha, ['fa_',num2str(namesGenerate(j))], TzerosArray);
    
    % Get all derivative
    T0d = [                    -cos(alpha(j)/2)*sin(alpha(j)/2) ,    exp(2*1i*p_rf)*sin(alpha(j)/2)*cos(alpha(j)/2) , -1i*exp(1i*p_rf) *cos(alpha(j)); ...
                exp(-2*1i*p_rf)*sin(alpha(j)/2)*cos(alpha(j)/2) ,                  -cos(alpha(j)/2)*sin(alpha(j)/2) ,  1i*exp(-1i*p_rf)*cos(alpha(j)); ...
            -1i/2*exp(-1i*p_rf)*cos(alpha(j))                   , 1i/2*exp(1i*p_rf)*cos(alpha(j))                   ,                  -sin(alpha(j))];
        
    % Get in cell arrays
    TdArray{j} = T0d;                                           % all derivatives     
    TdArray_partial_alpha.(['fa_',num2str(namesGenerate(j))]){j} = T0d;  % Derivative of T with respect to alpha(j)    
end   


% 1.3 - Organize
auxT_RFexc = -1i*exp(1i*p_ex)*sin(alpha_RF);

T      = blkdiag(auxT_RFexc, TArray{:}); % sin alfa
Td     = blkdiag(0, TdArray{:});

for j=1:n % loop over number of echoes  
    aux              = blkdiag(auxT_RFexc,TdArray_partial_alpha.(['fa_',num2str(namesGenerate(j))]){:});
    Td_partial_alpha = setfield(Td_partial_alpha, ['fa_',num2str(namesGenerate(j))], aux);
end
    

%% 2 - Loop over different relaxation rates

grad_a = [];        % Implementing the gradients according to NLO
phi_1_1 = [];       % dE / dalpha_l
phi_2_1 = [];       % d(dE_dT2) / dalpha_l
dPhi_by_z_1_1 = [];
dPhi_by_z_2_1 = [];


for iRate=1:nRates % loop over different relaxation rates - different T2 values
    
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

    % _---------------------- NLO -----------------_
    % NLO - 2.1 - derivatives of beta = TE = tau
    R0_beta  = diag([-R2*exp(-tau*R2)/2,-R2*exp(-tau*R2)/2,-R1*exp(-tau*R1)/2]);
    R0d_beta = diag(-exp(-tau*R2)/2+tau*R2*exp(-tau*R2)/2*[1 1 0]);  % derivative w.r.t R2 and beta    
    
    RArray_beta = cell(1,n);
    RArray(:)   = {sparse(R0_beta)};
    R_beta      = blkdiag(-R2*exp(-tau*R2)/2,RArray_beta{:});
    RArray(:)   = {sparse(R0d_beta)};
    Rd_beta     = blkdiag(-exp(-tau*R2)/2+tau*R2*exp(-tau*R2)/2,RArray{:});
    
    % derivatives of beta = TE = tau
    P_beta  = (R_beta*S);   % = RS
    Pd_beta = (Rd_beta*S);        
    % _--------------------------------------------_
    
    
    
    % Matrix representing the inter-echo duration
    E   = P*T*P;               % Transition Matrix no derivative
           
    % Transition Matrix with derivatives for T2
    EdR = Pd*T*P + P*T*Pd;     % Transition Matrix with derivatives for T2             
    
    
    
    % _---------------------- NLO -----------------_   
    % NLO - 2.2 - Get d (E) / d_partial_alpha (Ed_alpha) & d (dE/dT2) /d_partial_alpha (EdR_alpha)
    Ed_alpha  = [];
    EdR_alpha = [];
    
    for j=1:n 
        aux_Ed_alpha  = P*Td_partial_alpha.(['fa_',num2str(namesGenerate(j))])*P;                    % aux for: d (E     ) / d_alpha_l
        aux_EdR_alpha = Pd*Td_partial_alpha.(['fa_',num2str(namesGenerate(j))])*P + P*Td_partial_alpha.(['fa_',num2str(namesGenerate(j))])*Pd; % d (dE/dT2) / d_alpha_l
        Ed_alpha      = setfield(Ed_alpha, ['fa_',num2str(namesGenerate(j))],aux_Ed_alpha);         % d (E     ) / d_alpha_l
        EdR_alpha     = setfield(EdR_alpha, ['fa_',num2str(namesGenerate(j))],aux_EdR_alpha);        % d (dE/dT2) / d_alpha_l
    end
  
    % figures
    figure(88)
    for j=1:n
        subplot(2,3,j)
        dataforImage = Ed_alpha.(['fa_',num2str(namesGenerate(j))]);
        imagesc(abs(dataforImage))
        hold on
        caxis([0 0.6])        
    end
    figure(89)
    for j=1:n
        subplot(2,3,j)
        dataforImage = EdR_alpha.(['fa_',num2str(namesGenerate(j))]);
        imagesc(abs(dataforImage))
        hold on
        caxis([0 0.01])
    end
    
    
    % NLO - 2.3 - Get dE_partial_beta (Ed_alpha) & dE/dT2_partial_beta (EdR_alpha)
    Ed_beta  = P_beta*T*P  + P*T*P_beta;    % d (E     ) / d_beta
    EdR_beta = Pd_beta*T*P + Pd*T*P_beta + ...
               P_beta*T*Pd + P*T*Pd_beta;   % d (dE/dT2) / d_beta


    % NLO - 2.4 - Create Matrix Phi
    phi{1,1} = E;
    phi{1,2} = 0;
    phi{2,1} = EdR;
    phi{2,2} = E;
    
    % NLO - 2.5 - dPhi_dalpha
    dphi_dalpha{1,1} = Ed_alpha;
    dphi_dalpha{1,2} = 0;
    dphi_dalpha{2,1} = EdR_alpha;
    dphi_dalpha{2,2} = Ed_alpha;
    
    % NLO - 2.6 - dPhi_dbeta
    dphi_dbeta{1,1} = Ed_beta;
    dphi_dbeta{1,2} = 0;
    dphi_dbeta{2,1} = EdR_beta;
    dphi_dbeta{2,2} = Ed_beta;
    
    % _ -------------------------------------------_
    
    
    
    % Vector of all states
    x    = zeros(size(R,1),1);
    x(1) = 1;
    
    % === RF Excitation ===
    x = T*x;    
    
    
    
    % _------------- NLO --------------------_    
    % NLO - 2.7 - Initialize z - position x & position y - step 5-c)
    z{1,1} = x;
    z{2,1} = zeros(1+n*3,1);
         
    % _--------------------------------------_
    
    
    
    % === First echo ===
    % -> Deriviatives of state vector
    xdashR = EdR*x;
    
    % -> Update derivative value
    ds_dT2(1,iRate) = xdashR(1);

    % -> Calculate new state for next interation
    x            = E*x;
    epg(1,iRate) = x(1);    
    
    
    % _------------- NLO --------------------_

    
    % NLO - 2.8 - Position z(1) - position x & position y
    z{1,1+1} = phi{1,1}*z{1,1};
    z{2,1+1} = phi{2,1}*z{1,1} + phi{2,2}*z{2,1};
    
    % NLO - 2.9 - Initialize dz_dalpha = a & dz_dbeta = b - position x & position y
    dz_dalpha{1,1,1} = dphi_dalpha{1,1}.(['fa_',num2str(namesGenerate(1))])*z{1,1};    % Apenas conta a derivada do primeiro FA, pois k=1. (entradas de dz_dalpha (posicao no vetor, k, l) - l é a derivada de alpha
    dz_dalpha{2,1,1} = dphi_dalpha{2,1}.(['fa_',num2str(namesGenerate(1))])*z{1,1} + ...
                        dphi_dalpha{2,2}.(['fa_',num2str(namesGenerate(1))])*z{2,1};    
    
    dz_dbeta{1,1}  = dphi_dbeta{1,1}*z{1,1};
    dz_dbeta{2,1}  = dphi_dbeta{2,1}*z{1,1}   + dphi_dbeta{2,2}*z{2,1};    
    
    % _--------------------------------------_

    
    
    % === Subsequent echoes ===
    for i=2:n % loop over number of echoes
        
        % --- Get y variable ---            
        % -> Calculate derivatives using the product rule
        xdashR = EdR*x + E*xdashR; % Y_k - ( deriv. matrix transit. by states )  +  ( matrix transit. by (deriv. matrix transit. by previous states) )
        
        % -> a) Update derivative value
        ds_dT2(i,iRate) = xdashR(1);
        
        
        % _------------- NLO --------------------_
        % NLO - 2.10 - update z
        z{1,i+1} = phi{1,1}*z{1,(i+1)-1};
        z{2,i+1} = phi{2,1}*z{1,(i+1)-1} + phi{2,2}*z{2,(i+1)-1};    % (1+i)-1 , pq se soma 1 pelo RF de excitaçao, mas tira-se um pois é a posicao anterior

        % NLO - 2.11 - update dz_dalpha & dz_dbeta - position x & position y
        for j=1:i % loop over the variable 
             if j == i
                 dz_dalpha{1,i,j} = dphi_dalpha{1,1}.(['fa_',num2str(namesGenerate(j))]) * z{1,(i+1)-1};                             % (entradas de dz_dalpha (posicao no vetor, k, l) - l é a derivada de alpha
                 dz_dalpha{2,i,j} = dphi_dalpha{2,1}.(['fa_',num2str(namesGenerate(j))]) * z{1,(i+1)-1}  + ...
                                        dphi_dalpha{2,2}.(['fa_',num2str(namesGenerate(j))]) * z{2,(i+1)-1};
                 
             else
                 dz_dalpha{1,i,j} = phi{1,1} * dz_dalpha{1,i-1,j}  + ...
                                     dphi_dalpha{1,1}.(['fa_',num2str(namesGenerate(j))]) * z{1,(i+1)-1};                             % (entradas de dz_dalpha (posicao no vetor, k, l) - l é a derivada de alpha
                 dz_dalpha{2,i,j} = phi{2,1} * dz_dalpha{1,i-1,j}  +  phi{2,2} * dz_dalpha{2,i-1,j}  +  ...
                                     dphi_dalpha{2,1}.(['fa_',num2str(namesGenerate(j))]) * z{1,(i+1)-1}  +  ...
                                     dphi_dalpha{2,2}.(['fa_',num2str(namesGenerate(j))]) * z{2,(i+1)-1};                                                
             end               
        end         
        
        dz_dbeta{1,i} = phi{1,1}*dz_dbeta{1,i-1} + dphi_dbeta{1,1}*z{1,(i+1)-1};
        dz_dbeta{2,i} = phi{2,1}*dz_dbeta{1,i-1} + phi{2,2}*dz_dbeta{2,i-1} + ...
            dphi_dbeta{2,1}*z{1,(i+1)-1} + dphi_dbeta{2,2}*z{2,(i+1)-1};
        % _--------------------------------------_
    
    
        
        % -> b) Calculate new state for next interation
        x = E*x; % X_k
        
        % -> c) Get EPG        
        epg(i,iRate) = x(1);
        
    end
    
end

% figures
% check
figure(90)
a=1;
for i=1:n
    for j=1:n
        subplot(6,6,a)
        dataforImage =dz_dalpha{2,i,j};
        imagesc(abs(dataforImage))
        hold on
        caxis([0 0.004])
        a=a+1;
    end
end



%% Teste - get derivatives first point

sizZ     = size(z,2);
vector_a = zeros(sizZ-1,sizZ-1);

% derivadas do qual fazem depender o k ao longo dos echos
alpha_k1 = [dz_dalpha{2,1,1}];
alpha_k2 = [dz_dalpha{2,2,1} dz_dalpha{2,2,2} ];
alpha_k3 = [dz_dalpha{2,3,1} dz_dalpha{2,3,2} dz_dalpha{2,3,3} ];
alpha_k4 = [dz_dalpha{2,4,1} dz_dalpha{2,4,2} dz_dalpha{2,4,3} dz_dalpha{2,4,4} ];
alpha_k5 = [dz_dalpha{2,5,1} dz_dalpha{2,5,2} dz_dalpha{2,5,3} dz_dalpha{2,5,4} dz_dalpha{2,5,5} ];
alpha_k6 = [dz_dalpha{2,6,1} dz_dalpha{2,6,2} dz_dalpha{2,6,3} dz_dalpha{2,6,4} dz_dalpha{2,6,5} dz_dalpha{2,6,6}];

% derivadas do qual fazem depender o l ao longo dos echos
alpha_l1 = [dz_dalpha{2,1,1} dz_dalpha{2,2,1} dz_dalpha{2,3,1} dz_dalpha{2,4,1} dz_dalpha{2,5,1} dz_dalpha{2,6,1}];
alpha_l2 = [dz_dalpha{2,2,2} dz_dalpha{2,3,2} dz_dalpha{2,4,2} dz_dalpha{2,5,2} dz_dalpha{2,6,2}];
alpha_l3 = [dz_dalpha{2,3,3} dz_dalpha{2,4,3} dz_dalpha{2,5,3} dz_dalpha{2,6,3}];
alpha_l4 = [dz_dalpha{2,4,4} dz_dalpha{2,5,4} dz_dalpha{2,6,4}];
alpha_l5 = [dz_dalpha{2,5,5} dz_dalpha{2,6,5}];
alpha_l6 = [dz_dalpha{2,6,6}];

for jj=2:sizZ
    %x & y;
    x_signal(jj-1) = z{1,jj}(1);
    y(jj-1)        = z{2,jj}(1);
   
    % dy/d_partialAlpha = a's
    for ii=1:jj-1 % loop over the variable
        k = jj-1;
        l = ii;
        aux_a = dz_dalpha{2,jj-1,ii}(1);
        vector_a(jj-1,ii) = aux_a;
    end    
    % dy/d_beta = b
    b(jj-1) = dz_dbeta{2,jj-1}(1); % dy/dbeta
   
end

% get derivatives
x_signal = x_signal';
y        = y';
a        = vector_a;
b        = b';


%% Gradients Test



%% 3 - Gradients for optimization

% --------- test ------------------
for k = 1:sizZ-1    
    df_dalpha(k) = theta.theta_beta  * y' * a(:,k);
end
df_dbeta  = theta.dtheta_beta * y' * y   +  theta.theta_beta * y' * b;

% ---------------------------------

% % %     % _------------- NLO --------------------_
% % %     % NLO - 2.12 - Defining U matrix and dF_dalphaL & dF_dbeta
% % %     sizeN  = 1+n*3;
% % %     U{1,1} = zeros(sizeN,sizeN);
% % %     U{1,2} = eye(sizeN);    
% % %     
% % %     dF_dalphaL = [];
% % %     dF_dbeta = [];
% % %     
% % %     % NLO - 2.13 - Gradient of df dalphaL
% % %     for k=1:n
% % %         for ll=1:k
% % %             dF_dalphaL{1,k,ll} = theta.theta_beta*z{1,k+1}*U{1,1}*U{1,1}*dz_dalpha{1,k,ll}; % k+1 fruto de haver um pulso de RFexc na primeira posicao do vetor z
% % %             dF_dalphaL{2,k,ll} = theta.theta_beta*z{2,k+1}*U{1,2}*U{1,2}*dz_dalpha{2,k,ll};
% % %         end
% % %     end
% % %     
% % %     % NLO - 2.14 - Gradient of f dbeta
% % %     for k=1:n
% % %         for ll=1:k        
% % %             dF_dbeta{1,k,ll} = theta.dtheta_beta *   z{1,k+1}*U{1,1}*U{1,1}  + ...
% % %                  theta.theta_beta * z{1,k+1}*U{1,1}*U{1,1} * dz_dbeta{1,k,ll}; % k+1 fruto de haver um pulso de RFexc na primeira posicao do vetor z
% % %             dF_dbeta{2,k,ll} = theta.dtheta_beta *   z{2,k+1}*U{1,2}*U{1,2}  + ...
% % %                  theta.theta_beta * z{2,k+1}*U{1,2}*U{1,2} * dz_dbeta{2,k,ll};
% % %         end
% % %     end
% % %     % _--------------------------------------_
% % %     
% % %     
 
end