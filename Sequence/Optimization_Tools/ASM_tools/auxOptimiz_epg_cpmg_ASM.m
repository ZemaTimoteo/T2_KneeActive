% function [s,phasediag,P] = epg_cpmg(refoc_pulse,etl,T1,T2,esp)
%
%	EPG Simulation of CPMG sequence.  First flip angle
%	is 90 about y axis, and others by default are about
%	x-axis (make refoc_pulse complex to change that).
%
%	refoc_pulse = refocusing flip angle or list (radians)
%	etl = echo train length, if single refoc_pulse is to be repeated.
%	T1,T2,esp = relaxation times and echo spacing (arb. units).
%
%	Note that if refoc flip angle is constant, less than pi, and etl > 1 the
%	first refocusing flip angle is the average of 180 and the desired
%	refocusing flip angle, as proposed by Hennig.
%
%	All states are kept, for demo purposes, though this 
%	is not really necessary.

% function [s,phasediag,P] = auxOptimiz_epg_cpmg_ASM(exc_pulse,exc_phase,refoc_pulse,ETL,T1,T2,dTE,hennig,refoc_phase)
function [s,grad] = auxOptimiz_epg_cpmg_ASM(exc_pulse,exc_phase,refoc_pulse,ETL,T1,T2,dTE,hennig,refoc_phase)

%% ... 0 - Set parameters ...
% -- Default parameters:  ETL = length(refoc_pulse)
if (length(ETL)==0) ETL = length(refoc_pulse); end;

if (length(refoc_pulse)==1) && (ETL > 1) && (abs(refoc_pulse)<pi) && (hennig==1)
  % -- 1st flip reduced trick (Hennig)
  refoc_pulse(2)=refoc_pulse(1);
  refoc_pulse(1)=(pi*exp(i*angle(refoc_pulse(2)))+refoc_pulse(2))/2;
end;

% -- Allows for similar Refoc_pulse 
if (ETL > length(refoc_pulse)) refoc_pulse(end+1:ETL) = refoc_pulse(end); end; 


%% ... 1 - Initial conditions ...
% abs(refoc_pulse)*(180/pi) % degrees

P      = zeros(3,2*ETL);	% Allocate all known states, 2 per echo.
P(3,1) = 1;                 % Initial condition/equilibrium.

Pstore = zeros(4*ETL,ETL);	% Store F,F* states, with 2 gradients per echo
Zstore = zeros(2*ETL,ETL);	% Store Z states, with 2 gradients per echo
s_test(1) = P(1,1);

%% ... 2 - 90 excitation ...
% P = epg_rf(P,pi/2,pi/2);            % Do 90 tip.
P = epg_rf(P,exc_pulse,exc_phase);    % Do ex_refoc_pulse tip with phase phi (radians)
s_test(2) = P(1,1);
%% ... 3 - Refocus pulse ... 
s = zeros(1,ETL);		% Allocate signal vector to store.

for ech=1:ETL
  P = epg_grelax(P,T1,T2,dTE/2,1,0,1,1);            % Why dividing by 2?   % -- E - Relaxation Operator - Left crusher
  P = epg_rf(P,abs(refoc_pulse(ech)),angle(refoc_pulse(ech)));   % -- R - Rotation Operator - Refoc. RF
  P = epg_grelax(P,T1,T2,dTE/2,1,0,1,1);                         % -- E - Relaxation Operator - Right crusher

  s(ech)        = P(1,1);                   % Signal is F0 state.
  s_test(ech+2) = P(1,1);                   % Signal is F0 state.
  Pstore(2*ETL:4*ETL-1,ech) = P(2,:).';     % Put in negative states
  Pstore(1:2*ETL,ech) = flipud(P(1,:).');   % Put in positive, overwrite center.
  Zstore(:,ech) = P(3,:).';                 % Get Z

end

plotstate = cat(1,Pstore,Zstore);
% % dispim(plotstate);
% % xlabel('Echo'); ylabel('F(top) and Z(bottom) States');
phasediag = plotstate;	



%% ... 4 - Calculate Adjoint States ...
% % % % % 
% % % % % % parameters
% % % % % initial = 3;             % focus only on states from t=initial; 3 is first echo
% % % % % klim    = 25;            % order of maximum k, use 25
% % % % % Ny      = 4; Nx = Ny;    % TODO - resolution
% % % % % B1      = ones(Nx*Ny,1); % TODO (add to function inputs) - Necessary shape Number of points - number of channels [Ns, Nch]
% % % % % Nch     = 1;             % TODO number of B1 coils
% % % % % [Ns, ~] = size(B1);      % res(Nx*Ny), number of channels
% % % % % Nt      = ETL + 1;       % number of RF pulses (1x RF exc + ETL x RF refoc)
% % % % % np      = Nt/(Nch);      % number of states per B1 coil channel
% % % % % 
% % % % % klim = min([klim, round((np-1)/2)]);      % value for maximum number of allowed configuration states
% % % % % 
% % % % % kmax = 2*np - 1;                % up to just before the next RF pulse
% % % % % kmax = min([2*np-1, 2*klim-1]);
% % % % % N    = 3*(kmax-1)/2;            % number of states in total
% % % % % 
% % % % % c              = zeros(1,Nt+1);
% % % % % c(initial:end) = 1; 
% % % % % 
% % % % % % Calculate RS - Relaxation by Shift
% % % % % RS = RScalc(dTE,T1,T2,N);
% % % % % 
% % % % % % initializartion
% % % % % LL          = zeros(np+1,N*2,Ns);  % initialization of adjoint states per image space
% % % % % RRSS        = [RS, RS*0;RS*0, RS]; % Expand Matrix of RS
% % % % % RRSS(N+1,:) = -RRSS(N+1,:);
% % % % % 
% % % % % % alpha & phase
% % % % % alph        = repmat([exc_pulse,refoc_pulse],Nx*Ny,1);
% % % % % ph          = zeros(Nx*Ny,np);
% % % % % ph(:,2:end) = ph(:,2:end) + pi/2; % add CPMG phase
% % % % % 
% % % % % % Calculate the other states
% % % % % for ns = 1:Ns % loop over image space
% % % % %     L = zeros(np+1,N*2); % L - adjoint state variable
% % % % %     for jj=np:-1:2 % loop over number of different states (top-down)
% % % % % 
% % % % %         if jj == 2 % first adjoint state
% % % % %             A = Trot_fun_ASM(alph(ns,jj+1),ph(ns,jj+1)); % Rotation matrix direct definition
% % % % %             T = build_T_matrix_sub_ASM(A,3);
% % % % %             
% % % % %             % relaxation matrix
% % % % %             d          = zeros(6,6);
% % % % %             d(1)       = exp(-0.5*dTE/T2);
% % % % %             d(3+1,3+1) = exp(-0.5*dTE/T2);
% % % % %             
% % % % %             L(jj,[1:3 N+1:N+3]) = ...  % L*A*RS
% % % % %                 L( jj+1 , [1:3 N+1:N+3] ) * [real(T) -imag(T);imag(T) real(T)] * d ...
% % % % %                 + ...  % C(n+1) * f(n)
% % % % %                 [0, c(jj+1) * ( s_test(jj+1) ), 0, 0, c(jj+1) * ( s(jj-1) ), 0];
% % % % %         
% % % % %         elseif jj == np % last adjoint state        
% % % % %             temp      = zeros( 1 , 2*N );   
% % % % %             temp(2)   = c(jj+1)*( s_test(jj+1) );   % eq. 17 - lambda(N-1) - ultimo estado do EPG - ultimo estado do Target EPG
% % % % %             temp(N+2) = c(jj+1)*( s_test(jj+1) );   % TODO understand better why the semetry
% % % % %             
% % % % %             L(jj,:)   = temp;  % adjoint-state
% % % % %             
% % % % %         else % middle adjoint state                        
% % % % %             A = Trot_fun_ASM(alph(ns,jj+1),ph(ns,jj+1)); % Rotation matrix direct definition
% % % % %             
% % % % %             % Get: Lambda(n-1)
% % % % %             
% % % % %             %Lambda kron R (rotation Matrix = A)
% % % % %             %blockdiag_mult: given A square matrix, it returns   ' w =  kron(A,eye(N))*v  '   efficiently
% % % % %             temp = [ ...      % rotation Matrix - P, lambda - L =
% % % % %                       blockdiag_mult_ASM( real(A).' , L(jj+1,1:N).' )  +  blockdiag_mult_ASM( imag(A).' , L(jj+1,N+1:end).') ;...
% % % % %                     - blockdiag_mult_ASM( imag(A).' , L(jj+1,1:N).' )  +  blockdiag_mult_ASM( real(A).' , L(jj+1,N+1:end).') ...
% % % % %                    ].';
% % % % %             
% % % % %             % Output (Lamvda(n-1)*kron*R) * (E*S)
% % % % %             temp = temp * RRSS;
% % % % %             
% % % % %             % f(n+1)^H*C(n+1)
% % % % %             temp1      = zeros(1,2*N);
% % % % %             temp1(2)   = c(jj+1)*( s_test(jj+1) );
% % % % %             temp1(N+2) = c(jj+1)*( s_test(jj+1) );
% % % % %             
% % % % %             L(jj,:) = temp+temp1;
% % % % % 
% % % % %         end
% % % % %     end
% % % % %     LL(:,:,ns) = L; 
% % % % % end
% % % % % 
% % % % % %% ... 5 - Calculate Gradients ...
% % % % % GRAD = zeros(2*np*Nch,Ns);
% % % % % 
% % % % % % Acc B1 in xx and yy
% % % % % xeff = real(alph); % effective real
% % % % % yeff = imag(alph); % effective imag
% % % % % 
% % % % % for ns = 1:Ns % loop over space
% % % % %     
% % % % %     for jj = 2:np+1 % loop over number of different states
% % % % %         % Get derivatives in order to alpha and to phase
% % % % %         dT_daeff = dTrot_fun_da_ASM( alph(ns,jj-1) , ph(ns,jj-1) ); 
% % % % %         dT_dpeff = dTrot_fun_dp_ASM( alph(ns,jj-1) , ph(ns,jj-1) );
% % % % %         
% % % % %         % get relaxation matrix        
% % % % %         d          = zeros(N,N);
% % % % %         d(1,1)     = exp(-0.5*dTE/T2);     
% % % % % 
% % % % %         % get states 
% % % % %         temp1  = d * s_test(1,jj-1)';       % Relaxation by state
% % % % %         temp2  = RS * s_test(1,jj-1)';    % Relaxation by Shift by State
% % % % %         
% % % % %         % lambda * delta(P)/delta(alfa) * state
% % % % %         TEMP1a = LL(jj-1,:,ns) * [  blockdiag_mult_ASM( real(dT_daeff) , temp1(1:N) )  -  blockdiag_mult_ASM( imag(dT_daeff) , temp1(N+1:end) ) ;...
% % % % %                                     blockdiag_mult_ASM( imag(dT_daeff) , temp1(1:N) )  +  blockdiag_mult_ASM( real(dT_daeff) , temp1(N+1:end) )  ];
% % % % %         TEMP1p = LL(jj-1,:,ns) * [  blockdiag_mult_ASM( real(dT_dpeff) , temp1(1:N) )  -  blockdiag_mult_ASM( imag(dT_dpeff) , temp1(N+1:end) ) ;...
% % % % %                                     blockdiag_mult_ASM( imag(dT_dpeff) , temp1(1:N) )  +  blockdiag_mult_ASM( real(dT_dpeff) , temp1(N+1:end) )  ];
% % % % %         TEMP2a = LL(jj-1,:,ns) * [  blockdiag_mult_ASM( real(dT_daeff) , temp2(1:N) )  -  blockdiag_mult_ASM( imag(dT_daeff) , temp2(N+1:end) ) ;...
% % % % %                                     blockdiag_mult_ASM( imag(dT_daeff) , temp2(1:N) )  +  blockdiag_mult_ASM( real(dT_daeff) , temp2(N+1:end) )  ];
% % % % %         TEMP2p = LL(jj-1,:,ns) * [  blockdiag_mult_ASM( real(dT_dpeff) , temp2(1:N) )  -  blockdiag_mult_ASM( imag(dT_dpeff) , temp2(N+1:end) ) ;...
% % % % %                                     blockdiag_mult_ASM( imag(dT_dpeff) , temp2(1:N) )  +  blockdiag_mult_ASM( real(dT_dpeff) , temp2(N+1:end) )  ];
% % % % %         
% % % % %         % Derivative per axxes
% % % % %         dadx = xeff(ns,jj-1)  / alph(ns,jj-1)   * real( B1(ns) )  +  yeff(ns,jj-1) / alph(ns,jj-1)   * imag( B1(ns) );
% % % % %         dpdx = xeff(ns,jj-1)  / alph(ns,jj-1)^2 * imag( B1(ns) )  -  yeff(ns,jj-1) / alph(ns,jj-1)^2 * real( B1(ns) );
% % % % %         dady = -xeff(ns,jj-1) / alph(ns,jj-1)   * imag( B1(ns) )  +  yeff(ns,jj-1) / alph(ns,jj-1)   * real( B1(ns) );
% % % % %         dpdy = xeff(ns,jj-1)  / alph(ns,jj-1)^2 * real( B1(ns) )  +  yeff(ns,jj-1) / alph(ns,jj-1)^2 * imag( B1(ns) );
% % % % % 
% % % % %         % Get Gradient
% % % % %         if jj == 3 % first refocusing
% % % % %             GRAD( jj-1 , ns )        = dadx * TEMP1a  +  dpdx * TEMP1p;
% % % % %             GRAD( np*Nch+jj-1 , ns ) = dady * TEMP1a  +  dpdy * TEMP1p;
% % % % %         else % excitation & other refocusing
% % % % %             GRAD( jj-1 , ns )        = dadx * TEMP2a  +  dpdx * TEMP2p; 
% % % % %             GRAD( np*Nch+jj-1 , ns ) = dady * TEMP2a  +  dpdy * TEMP2p; 
% % % % %         end
% % % % %                     
% % % % %     end
% % % % %     
% % % % % end
% % % % % 
% % % % % grad = real(sum(GRAD,2));
% % % % % grad = grad(:);
% % % % % 
% % % % % %% 9 - with respect to parameters
% % % % % A = zeros(length(frequencies),sum(frequencies));
% % % % % for j = 2:length(frequencies)
% % % % %     A(j,sumfreq(j-1)+1:sumfreq(j-1)+frequencies(j)) = 1;
% % % % % end
% % % % % A(1,1:frequencies(1)) = 1;
% % % % % A = sparse(A);
% % % % % A = kron(eye(2*Nch),A);
% % % % % grad = A*grad;

end
