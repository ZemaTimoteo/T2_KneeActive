
function RS = RScalc(ESP,T1,T2,N)
%% Shift Matrix by Relaxation Matrix


%% 1 - Shift matrix ==== build Shift matrix, S with indices 
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
%S(end-1,end-2) = 1; % not do this, it influences the truncated version!
S = sparse(S);

%% 2 - Relaxation =====: note, Z0 regrowth not handled (No Z0 state here)
E1 = exp(-ESP/T1);
E2 = exp(-ESP/T2);
%R = diag(repmat([E2 E2 E1],[1 kmax+1]));
R                        = eye(N);
ii                       = 1:(N/3);
R(3*N*(ii-1)+3*(ii-1)+1) = E2;
R(3*N*ii-2*N+3*(ii-1)+2) = E2;
R(3*N*ii-N+3*(ii-1)+3)   = E1;

%% 3 - composites
RS=sparse(R*S);


