%% EPG forward for FSE

function [F, grad] = auxOptimiz_EPG_forward1(theta,Pha,varargin)

% F0 is the FID/Z0 created by the flipback pulse. Initialize in case not set
F0 = 0;

np   = length(theta);   % ETL
kmax = 2*np - 1;        % Number of total k (y-axis) - up to just before the next RF pulse
N    = 3*(kmax-1)/2;    % number of states in total - all kmax + the starting point for each pulse

% split magnitude and phase
alpha = abs(theta);
phi   = Pha;

% add CPMG phase
%alpha(2:end) = alpha(2:end)*exp(1i*pi/2);
% % phi(2:end) = phi(2:end) + pi/2;  % test 29/12/21

%% get variables
T1  = inf;
T2  = inf;
ESP = 10;
Npathway = inf;

klimit   = false;
flipback = false;

for ii=1:length(varargin)
    
    if strcmp(varargin{ii},'T1')
        T1 = varargin{ii+1};
    end
    if strcmp(varargin{ii},'T2')
        T2 = varargin{ii+1};
    end
    if strcmp(varargin{ii},'target')
        target = varargin{ii+1};
    end
    if strcmp(varargin{ii},'Fend')
        Fend = varargin{ii+1};
    end
    if strcmp(varargin{ii},'ESP')||strcmp(varargin{ii},'TE')
        ESP = varargin{ii+1};
    end
    % # of coherence pathways to consider (default is inf)
    if strcmp(varargin{ii},'Npathway')||strcmp(varargin{ii},'npath')
        Npathway = varargin{ii+1};
    end
    % more drastic version of above (see lab book 9-7-11)
    if strcmp(varargin{ii},'klimit')
        klimit=true;
    end
    % Allow user to define starting Equilibrium magnetization
    if strcmp(varargin{ii},'E')||strcmp(varargin{ii},'M0')
        E = varargin{ii+1};
    end
    % Simulate Flip-back pulse at the time of the last echo
    if strcmpi(varargin{ii},'flipback')||strcmpi(varargin{ii},'drive')
        flipback = true;
        alpha_fb = varargin{ii+1}; % COMPLEX
    end
end

% enforce pathway limit
if (N>Npathway)&&~klimit
    N=Npathway;
end

if klimit
    nr = length(theta)-1; % number of refocus pulses
    kmax = nr - 1 + mod(nr,2);
    N = 1.5*(nr+mod(nr,2)); % number of states in total
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
    NMAX = 3:3:3*(np-1);
    NMAX(NMAX>N)=(N-mod(N,3));
end
    
    
%% Shift ==========  S : build Shift matrix, S with indices 
S = zeros([N N]);
%%% F(k>1) look @ states just BEFORE np+1 pulse
kidx = 4:3:N; % miss out F1+
sidx = kidx-3;
%idx = sub2ind([N N],kidx,sidx);
idx = kidx + N*(sidx-1);
S(idx)=1;

%%% F(k<1) 
kidx = 2:3:N;
kidx(end)=[];% most negative state relates to nothing; related row is empty
sidx = kidx+3;
% ix = sub2ind([N N],kidx,sidx);
ix = kidx + N*(sidx-1);
S(ix)=1;

%%% Z states
kidx = 3:3:N;
%ix = sub2ind([N N],kidx,kidx);
ix = kidx + N*(kidx-1);
S(ix)=1;

%%% finally F1+ - relates to F-1-
S(1,2)=1;
S(end-1,end-2) = 1;

%% Relaxation =====  E : note, Z0 regrowth not handled (No Z0 state here) 
E1 = exp(-0.5*ESP/T1);  % dividing by to i.o.t. be similar of epg_cpmg_ASM - to split before and after RF 
E2 = exp(-0.5*ESP/T2);

%% Rotation =======  R
%R = diag(repmat([E2 E2 E1],[1 kmax+1]));
R = eye(N);
ii = 1:(N/3);
R(3*N*(ii-1)+3*(ii-1)+1)=E2;
R(3*N*ii-2*N+3*(ii-1)+2)=E2;
R(3*N*ii-N+3*(ii-1)+3)=E1;

%%%% composites
RS=R*S;

%% F matrix (many elements zero, not efficient)
F = (zeros([N np+1])); %% records state 

% initial state:
F(3,1) = 1; 

% Now the other states
for jj=2:np+1 
    % excitation
    if jj == 2 
        % Excitation pulse
        A         = Trot_fun_ASM(alpha(jj-1),phi(jj-1)); %<---- Rotation matrix direct definition
        F(1:3,jj) = A*F(1:3,jj-1);                       %<---- state straight after excitation [F0 F0* Z0]

    % refocusing
    else        
        % NMAX is maximum state index
        kidx = 1:NMAX(end);
        
        % Relaxation before RFocus
        temp       = RS*F(kidx,jj-1);   
        
        % Rotation with RFocus       
        A          = Trot_fun_ASM(alpha(jj-1),phi(jj-1));         % Rotation matrix direct definition
        temp(1)    = conj(temp(1));
        T          = build_T_matrix_sub_ASM(A,NMAX(end));
        F(kidx,jj) = T*temp;
        
        % Relaxation after RFocus         
        F(kidx,jj) = RS*F(kidx,jj);           
    end 
end

%% now adjoint states
if nargout > 1 % need gradient
    C = zeros(N,N);
    C(2,2) = 1; % sampling matrix
    C = sparse(C);
    L = sparse(zeros(np+1,N));
    
    % Now the other states
    for jj=np:-1:1
        if jj == 1 %
            A = Trot_fun(alpha(jj+1),phi(jj+1));    % Rotation matrix direct definition
            T = build_T_matrix_sub(A,3);            % reduced version
            d = zeros(3,3);
            d(1) = exp(-0.5*ESP/T2);d = sparse(d);
            L(jj,1:3) = L(jj+1,1:3)*T*d+[0,conj(F(2,jj+1)-target(2,jj+1)),0];
        elseif jj == np
            L(jj,:) = (F(:,jj+1)-target(:,jj+1))'*C;
        else
            
            A = Trot_fun(alpha(jj+1),phi(jj+1));
            % NMAX is maximum state index
            kidx = 1:NMAX(end);
            T = build_T_matrix_sub(A,NMAX(end));
            temp = L(jj+1,:)*T*RS;
            %temp(1) = conj(temp(1));
            % First evolve, then conj, then flip
            L(jj,:) = temp+(F(:,jj+1)-target(:,jj+1))'*C;
            
        end
    end
    
    %% now gradient
    grad = zeros(np,1);
    for jj = 2:np+1
        dT = build_T_matrix_sub(dTrot_fun_da(alpha(jj-1),phi(jj-1)),NMAX(end));
        if jj == 3
            d = zeros(N,N);
            d(1) = exp(-0.5*ESP/T2);d = sparse(d);
            grad(jj-1) = real(L(jj-1,:)*dT*d*F(:,jj-1));
        else
            grad(jj-1) = real(L(jj-1,:)*dT*RS*F(:,jj-1));
        end
    end
end

end

