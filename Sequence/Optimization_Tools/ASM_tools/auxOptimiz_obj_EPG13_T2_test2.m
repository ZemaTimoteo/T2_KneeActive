function [obj, grad, FF]  = auxOptimiz_obj_EPG13_T2_test2(params,ESP,T1,T2,Adj_stat,B1,frequencies,klim, testPlot)


%% EPG forward for FSE - objective function. maximizes the signal

% klim limit maximum k coefficient
%
% parameters are with respect to frequencies
%
% efficient implementation of blockdiagonal matrix vec multiplication
% w.r.t. obj_EPG8
%
% c:      is vector of {0,1} samplings, for instance: c(t) = 1 counts, c(t) = 0
%         does not count.
%         includes optimization w.r.t. phi (phase)
%         includes 2D spatially resolved , Nch channels
% B1 :    Nvoxels x Nch array
% params: [real(theta1) real(theta2) ... imag(theta1) imag(theta2) ...];
%         where index stands for channel
%
% NOTE: objective w.r.t. real and imaginary parts of theta!
%
% Adapted from Alessandro Sbrizzi March 2015 (based on Shaihan Malik EPG implementation)
% by: TTFernandes, IST - Dez, 2021

%% 1 - construct angles -> params2
[Ns, Nch] = size(B1); % res, number of channels
np        = length(params)/(2*Nch);
x         = reshape(params(1:np*Nch),np,Nch).'; %real parts
y         = reshape(params(np*Nch+1:end),np,Nch).'; % imaginary parts

params0   = x +1i*y;
params1   = zeros(Nch,sum(frequencies));
sumfreq   = cumsum(frequencies);

for j = (length(frequencies)):-1:2
    params1(:,sumfreq(j-1)+1:sumfreq(j)) = repmat(params0(:,j),[1 frequencies(j)]);
end

params1(:,1:frequencies(1)) = repmat(params0(:,1),[1 frequencies(1)]);
params1 = params1.';
params2 = [real(params1(:).') , imag(params1(:).')];
params  = params2;


%% 2 -
[Ns, Nch] = size(B1);
F0  = 0;
np  = length(params)/(2*Nch);

th  = zeros(Ns,np,Nch);
th  = B1*(x+1i*y);      % effective theta - multiple channels

alph     = [zeros(Ns,1) abs(th)];     % effective alpha
ph       = [zeros(Ns,1) angle(th)];   % effective phi
ph (:,2) = pi/2;

% add CPMG phase
xeff     = [zeros(Ns,1) real(th)];    % real theta


kmax = 2*np - 1;                    % up to just before the next RF pulse
kmax = min([2*np-1, 2*klim-1]);
N    = 3*(kmax-1)/2;                % number of states in total - F+; F-; Mz p/echo

% split magnitude and phase
Npathway = inf;
klimit   = false;

% enforce pathway limit
if (N>Npathway)&&~klimit
    N=Npathway;
end

if klimit
    nr   = size(th,2)-1; % number of refocus pulses
    kmax = nr - 1 + mod(nr,2);
    N    = 1.5*(nr+mod(nr,2)); % number of states in total
    if mod(nr,2) 
        % odd
        KMAX = [1:2:kmax (kmax-2):-2:1];
    else
        %even
        KMAX = [1:2:kmax kmax:-2:1];
    end
    NMAX = 1.5*(KMAX+1);
else
    % previous implementation
    NMAX         = 3:3:3*(np-1);
    NMAX(NMAX>N) = (N-mod(N,3));
end
    
    
%% 3 - Shift matrix ==== build Shift matrix, S with indices 
S = zeros([N N]);

%%% F(k>1) look @ states just BEFORE np+1 pulse
kidx   = 4:3:N; % miss out F1+
sidx   = kidx-3;
idx_N    = kidx + N*(sidx-1);
S(idx_N) =1;

%%% F(k<1) 
kidx      = 2:3:N;
kidx(end) = [];% most negative state relates to nothing; related row is empty
sidx      = kidx+3;
ix        = kidx + N*(sidx-1);
S(ix)     = 1;

%%% Z states
kidx  = 3:3:N;
ix    = kidx + N*(kidx-1);
S(ix) = 1;

%%% finally F1+ - relates to F-1-
S(1,2) = 1;
S      = sparse(S);

%% 4 - Relaxation =====: note, Z0 regrowth not handled (No Z0 state here)
E1 = exp(-0.5*ESP/T1); % dividing ESP by to i.o.t. be similar of epg_cpmg_ASM - to split before and after RF 
E2 = exp(-0.5*ESP/T2); % dividing ESP by to i.o.t. be similar of epg_cpmg_ASM - to split before and after RF 

R_matrix                        = eye(N);
ii                              = 1:(N/3);
R_matrix(3*N*(ii-1)+3*(ii-1)+1) = E2;
R_matrix(3*N*ii-2*N+3*(ii-1)+2) = E2;
R_matrix(3*N*ii-N+3*(ii-1)+3)   = E1;

% Full RS matrix
E1_full = exp(-ESP/T1); % dividing ESP by to i.o.t. be similar of epg_cpmg_ASM - to split before and after RF 
E2_full = exp(-ESP/T2); % dividing ESP by to i.o.t. be similar of epg_cpmg_ASM - to split before and after RF 

R_full                        = eye(N);
ii                            = 1:(N/3);
R_full(3*N*(ii-1)+3*(ii-1)+1) = E2_full;
R_full(3*N*ii-2*N+3*(ii-1)+2) = E2_full;
R_full(3*N*ii-N+3*(ii-1)+3)   = E1_full;


%% 5 - composites
RS      = R_matrix*S;
RS_full = R_full*S;

%% 6 - F matrix (many elements zero, not efficient)
FF = zeros(N*2,np+1,Ns);

% Now the other states
for ns = 1:Ns % loop over space
    % Initialization
    F = zeros([N*2 np+1]); %% records state 
    % initial state:
    F(3,1) = 1; 
    
    for jj = 2:np+1 %loop over time
        
        % excitation
        if jj == 2 
            % Excitation pulse
            A                   = Trot_fun_ASM(alph(ns,jj),ph(ns,jj)); %<---- Rotation matrix direct definition
            F([1:3 N+1:N+3],jj) = [real(A); imag(A)] * F(1:3,jj-1);    %<---- state straight after excitation [F0 F0* Z0]

        % refocusing            
        else
            % Relaxation before RFocus
            temp1    = RS * F(1:N,jj-1);
            temp2    = RS * F(N+1:2*N,jj-1);
            temp2(1) = -temp2(1);% D*temp;
            
            % Rotation with RFocus
            A = Trot_fun_ASM(alph(ns,jj),ph(ns,jj)); % Rotation matrix direct definition
            F(:,jj)  = [ blockdiag_mult_ASM(real(A),temp1) - blockdiag_mult_ASM(imag(A),temp2); ...
                         blockdiag_mult_ASM(imag(A),temp1) + blockdiag_mult_ASM(real(A),temp2) ];
            
            % Relaxation after RFocus
            F(1:N,jj)     = RS * F(1:N,jj);
            F(N+1:2*N,jj) = RS * F(N+1:2*N,jj);            

        end
    end
    FF(:,:,ns) = F; 
end

% Get Objective function
state = FF(2,:,:) + 1i*FF(N+2,:,:);
state = diag(Adj_stat) * squeeze(state);
obj   = 0.5*norm(state(:))^2;     % why? 29/12/21



% % figure()
% % plot(abs(FF(1,:,1)))
% % hold on, title('EPG States')



if nargout > 1 % need gradient

%% Get derivatives for dS/dT2

% Get derivatives of relaxation matrix in order to T2
TE   = ESP;
S0   = 1;     % Magnitude of B0

%%% ================== Math ==============================================
% S  = SUM( Fn * Cn )
% Fn = P(n-1) * F(n-1)
% P  = I (x) (RE)S

% dS/dT2 = SUM{n=1 to N-1}     [ PI{m=N-1 to n+1}[ P(m) ] ]    *    [ dP/dT2(n) ]   *   [F(n)]
%%% ======================================================================

% ================== get P(m) ===================================

for jj = 1:np+1 %loop over time    
    A      = Trot_fun_ASM(alph(1,jj),ph(1,jj)); % Rotation matrix direct definition
    if jj ==1      % Zero state
        P{jj} = zeros(3,3);
    elseif jj == 2 % Excitation %<--- P  = R
        P{jj} = A;
    else           % Refocusing %<--- P  = I (x) (RE)S
        P{jj}  = A*RS_full(1:3,1:3);
    end
end

% ================== get dP/dT2(n) =============================
% get derivative relaxation matrix w.r.t T2
dE_dT2eff      = dErelax_fun_dT2_ASM(T2,TE/2,S0,N); % derivative w.r.t T2
dRS            = dE_dT2eff*S;

dE_dT2eff_full = dErelax_fun_dT2_ASM(T2,TE,S0,N); % derivative w.r.t T2
dRS_full       = dE_dT2eff_full*S;

for jj = 1:np+1 %loop over time    
    A          = Trot_fun_ASM(alph(1,jj),ph(1,jj)); % Rotation matrix direct definition
    if jj ==1      % Zero state
        dP_dT2{jj} = zeros(3,3);
    elseif jj == 2 % derivative of Excitation %<--- dP/dT2  = dR/dT2 = 0
        dP_dT2{jj} =  zeros(3,3);
    else           % derivative of Refocusing %<--- dP/dT2  = I (x) (R * dE/dT2)S
        dP_dT2{jj} = A*dRS_full(1:3,1:3);
    end
end

clear dE_dT2eff dE_dT2eff_full


% ================== get dS/dT2(n) with all states ===========
%  get cycle of m variable 
%  get cycle of n variable 
% dS/dT2 = SUM{n=1 to N-1}     [ PI{m=N-1 to n+1}[ P(m) ] ]    *    [ dP/dT2(n) ]   *   [F(n)]

% initialize
nn = N;
Matrix_dS_dT2 = zeros(nn,np+1);

% for ns = 1:Ns % loop over space
    for N = 2:np+2 % loop over number of different states  
        idx_N = N-1;
                                     
        if N==2
            Matrix_dS_dT2(:,N) = blockdiag_mult_ASM(dP_dT2{N-1},FF(1:nn,N-1,ns)) ;
            
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
                aux_Matrix = blockdiag_mult_ASM(P_aux,blockdiag_mult_ASM(dP_dT2{n},FF(1:nn,n,ns)));
                aux_dS_dT2 = aux_Matrix + aux_dS_dT2 ;
                
                clear P_aux 
            end
            
            % Add P(m) independent part
            Matrix_dS_dT2(:,N) = blockdiag_mult_ASM(dP_dT2{N-1},FF(1:nn,N-1,ns)) + aux_dS_dT2 ;
            test_dS_dT2(:,N)   = aux_dS_dT2;
            clear aux_dS_dT2 
        end
    end
% end


% ================== get Gradient ===============================
grad = Matrix_dS_dT2(:,2:end);

% ================== plot =======================================
if testPlot == 'True'
    figure()
    plot(abs(Matrix_dS_dT2(2,size(Matrix_dS_dT2,2)-(size(alph,2)-2)+1:end)))
    hold on, title('Gradient')
    
end


end
end
