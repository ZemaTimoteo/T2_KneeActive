 function [obj, grad, FF]  = auxOptimiz_obj_EPG13_T2(params,ESP,T1,T2,c,B1,frequencies,klim,sigma)


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

alph     = abs(th);     % effective alpha
ph       = angle(th);   % effective phi
ph (:,1) = pi/2;

% add CPMG phase
xeff     = real(th);    % real theta


kmax = 2*np - 1;                    % up to just before the next RF pulse
kmax = min([2*np-1, 2*klim-1]);
N    = 3*(kmax-1)/2;                % number of states in total

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
idx    = kidx + N*(sidx-1);
S(idx) =1;

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

R                        = eye(N);
ii                       = 1:(N/3);
R(3*N*(ii-1)+3*(ii-1)+1) = E2;
R(3*N*ii-2*N+3*(ii-1)+2) = E2;
R(3*N*ii-N+3*(ii-1)+3)   = E1;

% Full RS matrix
E1_full = exp(-ESP/T1); % dividing ESP by to i.o.t. be similar of epg_cpmg_ASM - to split before and after RF 
E2_full = exp(-ESP/T2); % dividing ESP by to i.o.t. be similar of epg_cpmg_ASM - to split before and after RF 

R_full                        = eye(N);
ii                            = 1:(N/3);
R_full(3*N*(ii-1)+3*(ii-1)+1) = E2_full;
R_full(3*N*ii-2*N+3*(ii-1)+2) = E2_full;
R_full(3*N*ii-N+3*(ii-1)+3)   = E1_full;

%% 5 - composites
RS = sparse(R*S);
RS_full = sparse(R_full*S);

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
            A                   = Trot_fun_ASM(alph(ns,jj-1),ph(ns,jj-1)); %<---- Rotation matrix direct definition
            F([1:3 N+1:N+3],jj) = [real(A); imag(A)] * F(1:3,jj-1);        %<---- state straight after excitation [F0 F0* Z0]

        % refocusing            
        else
            % Relaxation before RFocus
            temp1    = RS * F(1:N,jj-1);
            temp2    = RS * F(N+1:2*N,jj-1);
            temp2(1) = -temp2(1);% D*temp;
            
            % Rotation with RFocus
            A = Trot_fun_ASM(alph(ns,jj-1),ph(ns,jj-1)); % Rotation matrix direct definition
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
state = diag(c) * squeeze(state);
obj   = 0.5*norm(state(:))^2;     % why? 29/12/21

% % figure()
% % plot(abs(FF(2,3:end,1)))

%% 7 - now adjoint states
if nargout > 1 % need gradient

% initializartion    
LL          = zeros(np+1,N*2,Ns);   % initialization of adjoint states per image space
RRSS        = [RS_full, RS_full*0;RS_full*0, RS_full];  % Expand Matrix of RS
RRSS(N+1,:) = -RRSS(N+1,:);

% Now the other states
for ns = 1:Ns            % loop over space
    L = zeros(np+1,N*2); % L - adjoint state variable
    
    for jj=np:-1:1       % loop over number of different states (top-down)

        if jj == 1 % first adjoint state
            A = Trot_fun_ASM(alph(ns,jj+1),ph(ns,jj+1)); % Rotation matrix direct definition
            T = build_T_matrix_sub_ASM(A,3);

            % relaxation matrix            
            d          = zeros(6,6);
            d(1)       = exp(-0.5*ESP/T2); % why? only F+, not F- neither Z?
            d(3+1,3+1) = exp(-0.5*ESP/T2);
            
            L(jj,[1:3 N+1:N+3]) = ... % L*A*RS
                L(jj+1,[1:3 N+1:N+3]) * [real(T) -imag(T);imag(T) real(T)] * d ...
                + ...
                [0, c(jj+1) * FF(2,jj+1,ns), 0, 0, c(jj+1) * FF(N+2,jj+1,ns), 0 ];
        
        elseif jj == np % last adjoint state        
            temp      = zeros(1,2*N);                 % N - number of states
            temp(2)   = c(jj+1) * FF(2,jj+1,ns);      % eq. 17 - lambda(N-1) - ultimo estado do EPG - 
            temp(N+2) = c(jj+1) * FF(N+2,jj+1,ns);    % TODO understand better why the semetry
            
            L(jj,:)   = temp;  % adjoint-state            
                        
        else % middle adjoint state

            A    = Trot_fun_ASM(alph(ns,jj+1),ph(ns,jj+1)); % Rotation matrix direct definition

            % Get: Lambda(n-1)
            
            %Lambda kron R (rotation Matrix = A)
            %blockdiag_mult: given A square matrix, it returns   ' w =  kron(A,eye(N))*v  '   efficiently
            temp = [ ...        % rotation Matrix - P, lambda - L =
                      blockdiag_mult_ASM( real(A).' , L(jj+1,1:N).' )  +  blockdiag_mult_ASM( imag(A).' , L(jj+1,N+1:end).') ;...
                    - blockdiag_mult_ASM( imag(A).' , L(jj+1,1:N).' )  +  blockdiag_mult_ASM( real(A).' , L(jj+1,N+1:end).') ...
                   ].';
            
            % Output (Lamvda(n-1)*kron*R) * (E*S)
            temp = temp * RRSS;

            % f(n+1)^H*C(n+1)   
            temp1      = zeros(1,2*N);
            temp1(2)   = c(jj+1) * FF(2,jj+1,ns);
            temp1(N+2) = c(jj+1) * FF(N+2,jj+1,ns);
            
            L(jj,:) = temp+temp1;

        end
    end
    LL(:,:,ns) = L; 
end

% % figure()
% % plot(abs(LL(3:end,2,1)))
% % imshow(LL(:,2,1),[])
%% 8 - now gradient

% Get derivatives of relaxation matrix in order to T2
TE   = [0 ESP/2 ESP:ESP:(np-1)*ESP];   % EchoTime for derivative
S0   = 1;     % Magnitude of B0
GRAD = zeros(2*np*Nch,Ns);
        
        
for ns = 1:Ns % loop over space
    
    for jj = 2:np+1 % loop over number of different states
        
        % get derivative relaxation matrix w.r.t T2
        dE_dT2eff    = dErelax_fun_dT2_ASM(T2,TE(jj-1),S0,N); % derivative w.r.t T2
        dRS_full     = sparse(dE_dT2eff*S);
        dRRSS        = [dRS_full, dRS_full*0;dRS_full*0, dRS_full];  % Expand Matrix of RS
        % dRRSS(N+1,:) = -dRRSS(N+1,:);
        if jj>2
            a = unique(dE_dT2eff);
            aux_dT2(jj-2) = a(2);
        end
        % get derivative relaxation matrix  for 1st RF pulse
        d          = zeros(2*N,2*N);
        d(1,1)     = S0*((0.5*ESP)/(T2.^2)).*exp(-0.5*ESP/T2);
        d(N+1,N+1) = S0*((0.5*ESP)/(T2.^2)).*exp(-0.5*ESP/T2);
        % get states
        temp1  = d     * FF(:,jj-1,ns);  % Relaxation by state
        temp2  = dRRSS * FF(:,jj-1,ns);  % Relaxation by Shift by State
        
        % get Rotation matrix direct definition
        A      = Trot_fun_ASM(alph(ns,jj-1),ph(ns,jj-1));

        % lambda * delta(P)/delta(alfa) * state        
        TEMP1a = LL(jj-1,:,ns) * [ blockdiag_mult_ASM(real(A),temp1(1:N)) - blockdiag_mult_ASM(imag(A),temp1(N+1:end));...
                                   blockdiag_mult_ASM(imag(A),temp1(1:N)) + blockdiag_mult_ASM(real(A),temp1(N+1:end)) ];
        TEMP1p = LL(jj-1,:,ns) * [ blockdiag_mult_ASM(real(A),temp1(1:N)) - blockdiag_mult_ASM(imag(A),temp1(N+1:end));...
                                   blockdiag_mult_ASM(imag(A),temp1(1:N)) + blockdiag_mult_ASM(real(A),temp1(N+1:end)) ];
        TEMP2a = LL(jj-1,:,ns) * [ blockdiag_mult_ASM(real(A),temp2(1:N)) - blockdiag_mult_ASM(imag(A),temp2(N+1:end));...
                                   blockdiag_mult_ASM(imag(A),temp2(1:N)) + blockdiag_mult_ASM(real(A),temp2(N+1:end)) ];
        TEMP2p = LL(jj-1,:,ns) * [ blockdiag_mult_ASM(real(A),temp2(1:N)) - blockdiag_mult_ASM(imag(A),temp2(N+1:end));...
                                   blockdiag_mult_ASM(imag(A),temp2(1:N)) + blockdiag_mult_ASM(real(A),temp2(N+1:end)) ];
        
        % Derivative per axxes            
        dadxj = xeff(ns,jj-1)  / alph(ns,jj-1)   * real(B1(ns)) ;
        dpdyj = xeff(ns,jj-1)  / alph(ns,jj-1)^2 * real(B1(ns)) ;

        % Get Gradient
        if jj == 3 % first refocusing
            GRAD(jj-1,ns)        = dadxj*TEMP1a;
            GRAD(np*Nch+jj-1,ns) = dpdyj*TEMP1p;
        else       % excitation & other refocusing
            GRAD(jj-1,ns)        = dadxj*TEMP2a; 
            GRAD(np*Nch+jj-1,ns) = dpdyj*TEMP2p; 
        end

    end
end

grad = real(GRAD);
grad = grad(:);

% % figure();plot(aux_dT2)
end
end
