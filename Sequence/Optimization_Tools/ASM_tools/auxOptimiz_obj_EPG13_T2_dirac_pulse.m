function [obj, grad, FF]  = auxOptimiz_obj_EPG13_T2_dirac_pulse(params,exc_puls,refoc_puls,ESP,ETL,T1,T2,Adj_stat,B1,frequencies,klim,Trec,Trec_Contrl,testPlot)


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

%% 6 - Run throw different moments of RF pulses

% 6.1 - ===== Initialize =====
% 6.1.1 - Initialization for F-states
FF        = zeros(N*2,np+1);
exc_pulse = size(params,1);
alph_test = [zeros(length(exc_puls),1) exc_puls refoc_puls];
ph_test   = [ph(1,:)];

% 6.1.2 - Get derivatives of relaxation matrix in order to T2
TE   = ESP;
S0   = 1;     % Magnitude of B0
alph_test = [zeros(length(exc_puls),1) exc_puls refoc_puls];
ph_test   = [ph(1,:)];

% ==== Cicle for points per RF pulse ----
for z=1:size(alph_test,1) % |-->pontos dos pulsos (128)
    % 6.2 - ==== F matrix (many elements zero, not efficient) ====
    % 6.2.1 - Get States
    s = auxOptimiz_epg_cpmg(np, alph_test(z,:), ph_test, RS, N);
    aux_soma(:,:,z) = s;

    % 6.3 - ====== Get derivatives for dS/dT2 ==========
    if nargout > 1 % need gradient
        % 6.3.1 - get dS/dT2(n) with all states
        gR = auxOptimiz_dEPG(alph_test(z,:), ph_test, RS_full, T2, TE, S0, N, S, np, s);
%         gR = auxOptimiz_dEPG_test(alph_test(z,:), ph_test, RS, T2, TE/2, S0, N, S, np, s);
        soma_dEPG(:,:,z) = gR;
    end
end

% 6.3 - Sum across all RF pulse coordenates (128 coordenates z)
F             = sum(aux_soma,3);
Matrix_dS_dT2 = sum(soma_dEPG,3);    % -- Somatorio ao longo da fatia, soma das 128 coordenadas z

%  figure; plot(abs(aux_soma(:,:)))
%  figure; plot(abs(soma_dEPG(:,:)))
clear aux_soma soma_dEPG

% 6.4 - normalize dictionary
% % FF = F ./ norm( F );        
FF = F;

% 6.5 - Get Objective function
state = FF(2,:) + 1i*FF(N+2,:);
state = diag(Adj_stat) * squeeze(state)';
obj   = 0.5*norm(state(:))^2;     % why? 29/12/21

% 6.6 - Get Gradient 
grad = Matrix_dS_dT2(:,2:end);

% 6.7 - Get Figure
if testPlot=='True'
    % 6.7.1 - Plot for Dictionary (F0 states)    
    figure()
    plot(abs(FF(2,3:end)))
    hold on, title('test EPG States')    
    % 6.7.2 - Plot for Gradient
    figure()
%     plot(abs(Matrix_dS_dT2(1,size(Matrix_dS_dT2,2)-(size(alph,2)-2)+1:end)))
    plot(abs(Matrix_dS_dT2(1,size(Matrix_dS_dT2,2)-(size(alph,2)-2)+1:end)))
    hold on, title('Gradient')    
end

%% 7 - Control for Trecovery & TR variation
if Trec_Contrl == 'True'
    auxVariab_FF     = exp( - Trec/T1);
    for ll=1:size(FF,1)
        aux2var_S_T2(ll) = (1 - auxVariab_FF) / ( 1 - auxVariab_FF * FF(ll,ETL+2) );
        FF(ll,:) = aux2var_S_T2(ll).* FF(ll,:);
    end
end



end
