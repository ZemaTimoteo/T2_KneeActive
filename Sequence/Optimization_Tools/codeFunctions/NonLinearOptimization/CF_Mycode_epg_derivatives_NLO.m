function [ grad ] = CF_Mycode_epg_derivatives_NLO( n, beta, T1, T2, alpha_RF, alpha, params, gamma)
                                
%% 1 - Define Variables and matrices
plotTest = params.plotTest;
tau      = beta/2;         % Half of echospacing (ESP/2) in (s) = TE/2

% ----- 1.1 - Get Rotation matrices (T) and dT_partial_alpha -----
p_rf = 0;
p_ex = pi/2;  

% 1.1.1 - RF excitation matrix, T (Rotation Matrix)
TArray{1}   =  [ (cos(alpha_RF/2))^2                  ,  exp(2*1i*p_ex)*(sin(alpha_RF/2))^2 , -1i*exp(1i*p_ex)*sin(alpha_RF);  ...
                  exp(-2*1i*p_ex)*(sin(alpha_RF/2))^2 , (cos(alpha_RF/2))^2                 ,  1i*exp(-1i*p_ex)*sin(alpha_RF); ...
                  -1i/2*exp(-1i*p_ex)*sin(alpha_RF)   ,  1i/2*exp(1i*p_ex)*sin(alpha_RF)    ,  cos(alpha_RF)];

% 1.1.2 - RF mixing matrix, T (Rotation Matrix)
for j=1:n % loop over number of echoes 
    Tzeros0 = zeros(3,3);
    T0   =  [ (cos(alpha(j)/2))^2                  ,  exp(2*1i*p_rf)*(sin(alpha(j)/2))^2 , -1i*exp(1i*p_rf)*sin(alpha(j));  ...
               exp(-2*1i*p_rf)*(sin(alpha(j)/2))^2 , (cos(alpha(j)/2))^2                 ,  1i*exp(-1i*p_rf)*sin(alpha(j)); ...
              -1i/2*exp(-1i*p_rf)*sin(alpha(j))    ,  1i/2*exp(1i*p_rf)*sin(alpha(j))    ,  cos(alpha(j))];
    % Ge in cell arrays   
    TArray{j+1}      = T0;
    TzerosArray{j} = Tzeros0;
end

% 1.1.3 - Derivative of T with respect to alpha
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


% ----- 1.2 - Define Functions of evolution -----
% Define E - Relaxation and shift Matrix and rotation
E      = @(aux_jj,aux_x,aux_T1,aux_T2,aux_tau)  shiftMatrix_NLO(  RelaxationMatrix_NLO(  TArray{aux_jj} * shiftMatrix_NLO(  RelaxationMatrix_NLO(aux_x,aux_T2,aux_T1,aux_tau))  ,aux_T2,aux_T1,aux_tau  )  );
% % E_test = @(aux_jj,aux_T1,aux_T2,aux_tau)  shiftMatrix_NLO(  RelaxationMatrix_NLO(  TArray{aux_jj} * shiftMatrix_NLO(  RelaxationMatrix_test_NLO(aux_T2,aux_T1,aux_tau))  ,aux_T2,aux_T1,aux_tau  )  );
 
% Define dE_dT2 
% RdS*T*RS + RS*T*RdS
dE_dT2 = @(aux_jj,aux_x,aux_T1,aux_T2,aux_tau)  shiftMatrix_NLO(  RelaxationMatrix_NLO(  TArray{aux_jj} * shiftMatrix_NLO(  RelaxationMatrix_dT2_NLO(aux_x,aux_T2,aux_T1,aux_tau))  ,aux_T2,aux_T1,aux_tau  )  ) ...
              + ...
              shiftMatrix_NLO(  RelaxationMatrix_dT2_NLO(  TArray{aux_jj} * shiftMatrix_NLO(  RelaxationMatrix_NLO(aux_x,aux_T2,aux_T1,aux_tau))  ,aux_T2,aux_T1,aux_tau  )  );
          
% ... 1.3 - Get d (E) / d_partial_alpha (Ed_alpha) & d (dE/dT2) /d_partial_alpha (EdR_alpha)
% Define E derivate to alpha partial - Pd*Td_partial_alpha*P + P*Td_partial_alpha*Pd                                                                                                                                                                                                    d
Ed_alpha  = @(aux_jj,aux_ii,aux_x,aux_T1,aux_T2,aux_tau)  shiftMatrix_NLO(  RelaxationMatrix_NLO(  TdArray_partial_alpha.(['fa_',num2str(namesGenerate(aux_jj-1))]){1,aux_ii-1} * shiftMatrix_NLO(  RelaxationMatrix_NLO(aux_x,aux_T2,aux_T1,aux_tau))  ,aux_T2,aux_T1,aux_tau  )  );
                                                                                                                                                                                                                                                                                    
% Define EdTE derivate to alpha partial - Pd*Td_partial_alpha*P + P*Td_partial_alpha.*Pd
EdR_alpha = @(aux_jj,aux_ii,aux_x,aux_T1,aux_T2,aux_tau)  shiftMatrix_NLO(  RelaxationMatrix_NLO(  TdArray_partial_alpha.(['fa_',num2str(namesGenerate(aux_jj-1))]){1,aux_ii-1} * shiftMatrix_NLO(  RelaxationMatrix_dT2_NLO(aux_x,aux_T2,aux_T1,aux_tau))  ,aux_T2,aux_T1,aux_tau  )  ) ...
              + ...                                                                                                                                                                                                                                                             
              shiftMatrix_NLO(  RelaxationMatrix_dT2_NLO( TdArray_partial_alpha.(['fa_',num2str(namesGenerate(aux_jj-1))]){1,aux_ii-1} * shiftMatrix_NLO(  RelaxationMatrix_NLO(aux_x,aux_T2,aux_T1,aux_tau))  ,aux_T2,aux_T1,aux_tau  )  );

% ... 1.4 - Derivatives in order to Beta (Beta/2 = Tau = TE/2)
% Derivative of P*T*P in order to beta  (Pbeta*T*P + P*T*Pbeta)
Ed_beta = @(aux_jj,aux_x,aux_T1,aux_T2,aux_tau)  ...
              shiftMatrix_NLO(  RelaxationMatrix_NLO(  TArray{aux_jj} * shiftMatrix_NLO(  RelaxationMatrix_dBeta_NLO(aux_x,aux_T2,aux_T1,aux_tau))  ,aux_T2,aux_T1,aux_tau  )  )  ...
              + ...
              shiftMatrix_NLO(  RelaxationMatrix_dBeta_NLO(  TArray{aux_jj} * shiftMatrix_NLO(  RelaxationMatrix_NLO(aux_x,aux_T2,aux_T1,aux_tau))  ,aux_T2,aux_T1,aux_tau  )  );

% Derivative of beta of dE_dT2 (Pd_beta*T*P + Pd*T*P_beta + P_beta*T*Pd + P*T*Pd_beta)
EdR_beta = @(aux_jj,aux_x,aux_T1,aux_T2,aux_tau)  ...
               shiftMatrix_NLO(  RelaxationMatrix_NLO(  TArray{aux_jj} * shiftMatrix_NLO(  RelaxationMatrix_dbeta_dT2_NLO(aux_x,aux_T2,aux_T1,aux_tau))  ,aux_T2,aux_T1,aux_tau  )  ) ...
               + ...
               shiftMatrix_NLO(  RelaxationMatrix_dBeta_NLO(  TArray{aux_jj} * shiftMatrix_NLO(  RelaxationMatrix_dT2_NLO(aux_x,aux_T2,aux_T1,aux_tau))  ,aux_T2,aux_T1,aux_tau  )  ) ...
               + ...
               shiftMatrix_NLO(  RelaxationMatrix_dT2_NLO(  TArray{aux_jj} * shiftMatrix_NLO(  RelaxationMatrix_dBeta_NLO(aux_x,aux_T2,aux_T1,aux_tau))  ,aux_T2,aux_T1,aux_tau  )  ) ...
               + ...
               shiftMatrix_NLO(  RelaxationMatrix_dbeta_dT2_NLO(  TArray{aux_jj} * shiftMatrix_NLO(  RelaxationMatrix_NLO(aux_x,aux_T2,aux_T1,aux_tau))  ,aux_T2,aux_T1,aux_tau  )  ) ;
          

%% 2 - Evolution of EPG states

% 2.1 - Initialization of states
clear z epg ds_dT2 dz_dbeta
x      = zeros(3,3);
y      = zeros(3,3);
x(3,1) = 1;                 % Initial condition/equilibrium.
ini_x  = ones(3,3);
ini_x(3,1) = 0;

% 2.2 - RF excitation 90º
x         = TArray{1}*x;
epg(1)    = x(1,1);
z{1}      = [x; y];


% 2.3 - Echoes (Relaxation and Shifting +  Rotation + Relaxation and Shifting)
for jj = 2:n+1
    % 2.3.1 - Get dz_dalpha_partial - position x & position y
    for ii=2:jj % loop over the variable
        if ii == jj
            dz_dalpha{1,jj,ii} = Ed_alpha(jj,ii,x,T1,T2,tau);                             % (entradas de dz_dalpha (posicao no vetor, k, l) - l é a derivada de alpha
            dz_dalpha{2,jj,ii} = EdR_alpha(jj,ii,x,T1,T2,tau)  + Ed_alpha(jj,ii,y,T1,T2,tau);
            
        else
            dz_dalpha{1,jj,ii} = E(jj,dz_dalpha{1,jj-1,ii},T1,T2,tau) + ...
                                    Ed_alpha(jj,ii,x,T1,T2,tau);     % (entradas de dz_dalpha (posicao no vetor, k, l) - l é a derivada de alpha
            dz_dalpha{2,jj,ii} = dE_dT2(jj,dz_dalpha{1,jj-1,ii},T1,T2,tau)   +   E(jj,dz_dalpha{2,jj-1,ii},T1,T2,tau)  ...
                                    +  ...
                                 EdR_alpha(jj,ii,x,T1,T2,tau)                +   Ed_alpha(jj,ii,y,T1,T2,tau);
        end
    end        
    
    % 2.3.2 - Get dz_dbeta
    dz_dbeta{1,jj} = Ed_beta(jj,x,T1,T2,tau);
    dz_dbeta{2,jj} = EdR_beta(jj,x,T1,T2,tau) + Ed_beta(jj,y,T1,T2,tau);

    % 2.3.3 - Get Updated y = EdR*x + E*y
    y   =  dE_dT2(jj,x,T1,T2,tau)+E(jj,y,T1,T2,tau);
    
    % 2.3.4 - Get Updated x
    x   =  E(jj,x,T1,T2,tau);
    
    % 2.3.5 - Get Updated z
    z{jj}  = [x; y];
    
    % 2.3.6 - get epg and ds_dT2
    epg(jj)      = x(1,1);
    ds_dT2(jj-1) = y(1,1);

end

%% figure
if plotTest == 'True'
    figure(90)
    a=1;
    for i=1:n
        for j=1:n
            subplot(6,6,a)
            dataforImage =dz_dalpha{2,i,j};
            imagesc(abs(dataforImage))
            hold on
            caxis([0 2])
            a=a+1;
        end
    end
end

%% 3 - Organize data
clear x y
sizZ     = size(z,2);
vector_a = zeros(sizZ-1,sizZ-1);

for jj=2:sizZ
    %x & y;
    x(jj-1) = z{1,jj}(1);
    y(jj-1) = z{1,jj}(4);
   
    % dy/d_partialAlpha = a's
    for ii=2:jj % loop over the variable
        k = jj; l = ii;
        aux_a = dz_dalpha{2,k,l}(1);
        vector_a(k-1,l-1) = aux_a;
    end    
    % dy/d_beta = b
    b(jj-1) = dz_dbeta{2,jj}(1); % dy/dbeta
   
end

% 3.1 - Get variables
x = x';
y = y';
a = vector_a;
b = b';

%% 4 - Gradients for optimization
% 4.1 - Grad ds_dT2
grad.ds_dT2 = real(ds_dT2);

% 4.2 - Grad dF_dalpha
for k = 1:sizZ-1    
    grad.dF_dalpha(1,k)   = real(gamma.gamma_beta  * y' * a(:,k) +  gamma.gamma_beta * a(:,k)' * y); % partial derivative for influence off all 'l=FA' 
% %     l=k;
% %     grad.aux_dF_dalpha(l) = real(gamma.gamma_beta  * y' * a(l,:)' +  gamma.gamma_beta * a(l,:) * y);
end

% 4.3 - Grad dF_dbeta
grad.dF_dbeta  = real(gamma.dgamma_beta * y' * y   + gamma.gamma_beta * b' * y + gamma.gamma_beta * y' * b);


