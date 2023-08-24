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

function [s_Fp,s_Fn,s_Z,phasediag,P] = epg_cpmg_1LA(exc_pulse,exc_phase,refoc_pulse,ETL,T1,T2,dTE,hennig,refoc_phase)

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

%% ... 2 - 90 excitation ...
% P = epg_rf(P,pi/2,pi/2);            % Do 90 tip.
P = epg_rf(P,exc_pulse,exc_phase);    % Do ex_refoc_pulse tip with phase phi (radians)

%% ... 3 - Refocus pulse ...
s = zeros(1,ETL);		% Allocate signal vector to store.

for ech=1:ETL
  P = epg_grelax(P,T1,T2,dTE/2,1,0,1,1);                         % -- Half TE relaxation before RF
% %   P = epg_grelax_1LA_noRelax(P,1,0,1,1);                         % -- Half TE relaxation before RF
  P = epg_rf(P,abs(refoc_pulse(ech)),angle(refoc_pulse(ech)));   % -- Refoc. RF
  P = epg_grelax(P,T1,T2,dTE/2,1,0,1,1);                         % -- Half TE relaxation after RF
% %   auxP = epg_grelax_1LA_noRelax(P,1,0,1,1);                         % -- Half TE relaxation after RF
% %   P = auxP(:,1:ETL*2);
  
  s_Fp(ech) = P(1,1);                          % Signal is F0 state.
  s_Fn(ech) = P(2,3);                          % Signal is F0 state.
% %   s_Fn(ech) = P(2,1);                          % Signal is F0 state.
  s_Z(ech)  = P(3,2);                          % Signal is F0 state.
% %   s_Z(ech)  = P(3,1);                          % Signal is F0 state.
  
  Pstore(2*ETL:4*ETL-1,ech) = P(2,:).';     % Put in negative states
  Pstore(1:2*ETL,ech) = flipud(P(1,:).');   % Put in positive, overwrite center.
  Zstore(:,ech) = P(3,:).';                 % Get Z
   



end
% % 
% % figure()
% % plot(abs(s))

plotstate = cat(1,Pstore,Zstore);
% % dispim(plotstate);
% % xlabel('Echo'); ylabel('F(top) and Z(bottom) States');
phasediag = plotstate;	
