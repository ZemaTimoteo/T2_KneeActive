 function [dictionary,pars] = dict_pars_generator(T1,T2,B1,ESP, ETL,refoc_phase,phase_exc,FA_exc,FA_refoc,dir_data)
% generation of a dictionary and its indexes (pars) with 3 entries: T2, the
% phase of the 90º RF (phi) and the values of B1.
% in each entry: the differents states of the magnetizations are present

%  Functions used:
% slr_slice
% my epg to generate the states over the slice excited (sum) for different
% T2
% epg.cmg to generate the states

% inputs:
% T1: constant- relaxation time of the longitudinal magnetization
% T2 : interval- relaxation time of the transversal magnetization
% B1_scale: interval- ???
% ESP: spacing time between echoes
% ETL: refocusing echoes train length
% refoc_phase = exp(zeros(1,ETL)./(180/pi).*1i); % phase of refoc pulses (0 = about x axis)
% phi = pi/2; % phase of excitation
% FA_exc - Flip angle excitation pulse: vector
% FA_refoc - Flip angle refocusing pulse: vector
% dTE - Spacing Time between Echoes (ms)

%% 1 - initialization of variable
exc_pulse       = [];
FA_refocUNIQUE  = unique(FA_refoc,'stable');
LargePhant_dict = zeros(ETL,length(T2)*length(phase_exc)*length(B1));

%% 2 - 'OR' slice profile with SLR method + Dictionary
    
% Loop over B1 scale values
for jj=1:size(B1,2)
% %     f = waitbar(jj/size(B1,2),'Building EPGs...');
    refoc_pulse     = [];
    
    % slice profile with SLR method
    for ii=1:size(FA_refocUNIQUE,2)
        [exc_pulse, refoc_pulseSLR] = slr_profile(B1(jj),FA_refocUNIQUE(ii),ESP,'HSM2',dir_data);   % corrige com os perfis adaptar para os pulsos da Siemens/pypulseq
        refoc_pulse = [refoc_pulse refoc_pulseSLR(:)];
        clear refoc_pulseSLR
    end
    
    Dict = my_epg(exc_pulse,refoc_pulse,phase_exc,refoc_phase,T1,T2,ESP,ETL);
    % add the entry to LargePhant_dic
    LargePhant_dict(:,(jj-1)*length(T2)*length(phase_exc)+1:jj*length(T2)*length(phase_exc)) = Dict;
    clear refoc_pulse

% %     close(f)
end


%% 3 - PARS generation
% Variable pars contains the values of T2, the phase of the 90º RF (phi)
% and the values of B1.

%initialization of variables:
pars_generation = zeros(length(T2)*length(phase_exc)*length(B1),3);

for j=1:length(T2)*length(phase_exc):size(pars_generation,1)
    ind_B1 = (j+length(T2)*length(phase_exc)-1)/(length(T2)*length(phase_exc));
    
    for i=j:length(T2):j+(length(T2)*length(phase_exc))-1
        pars_generation(i:i+length(T2)-1,1) = T2';
        ind_phi                             = (i+length(T2)-1)/length(T2)-((ind_B1-1)*length(phase_exc));
        pars_generation(i:i+length(T2)-1,2) = phase_exc(ind_phi)*ones(length(T2),1);
    end
    
    pars_generation(j:j+(length(T2)*length(phase_exc))-1,3) = B1(ind_B1)*ones(length(T2)*length(phase_exc),1);
end


%% 4 - Out parameters

dictionary = LargePhant_dict ;
pars = pars_generation ;

end

