function [dictionary,pars] = dict_pars_generator_EPGtest_Simulation(T1,T2,B1,ESP, ETL,refoc_phase,phase_exc,FA_exc,FA_refoc, params, methodDic)
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
if methodDic == 'DICT_SLR'
    
    % Loop over B1 scale values
    for jj=1:size(B1,2)
        f = waitbar(jj/size(B1,2),'Building EPGs...');
        refoc_pulse = [];
        rf_exc = FA_exc;
        
        % slice profile with SLR method
        for ii=1:size(FA_refocUNIQUE,2)
            [exc_pulse, refoc_pulseSLR] = slr_profile_wRFdesign(rf_exc,B1(jj),FA_refocUNIQUE(ii),params)
            refoc_pulse = [refoc_pulse refoc_pulseSLR(:)];
        end    
        
% %         %
% %         figure();plot(exc_pulse);hold on; yline(max(exc_pulse)/2)
% %         figure();plot(refoc_pulse);hold on; yline(max(refoc_pulse)/2)        
% %         %
        
        Dict = my_epg(exc_pulse,refoc_pulse,phase_exc,refoc_phase,T1,T2,ESP,ETL);
        % add the entry to LargePhant_dic
        LargePhant_dict(:,(jj-1)*length(T2)*length(phase_exc)+1:jj*length(T2)*length(phase_exc)) = Dict;
        
        clear refoc_pulse
    end
    
%% 2 - 'OR' Just generate the dictionary
elseif methodDic == 'JUST_DIC'
    % initialization of variable
    LargePhant_dict = zeros(ETL,length(T2)*length(phase_exc)*length(B1));
    
    % Loop over B1 scale values
    for i=1:length(B1)
        f = waitbar(i/size(B1,2),'Building EPGs...');

        exc_pulse(i) = FA_exc*B1(i);
        refoc_pulse  = FA_refoc*B1(i);
        Dict = my_epg(exc_pulse(i),refoc_pulse,phase_exc,refoc_phase,T1,T2,ESP,ETL);
        % add the entry to LargePhant_dic
        LargePhant_dict(:,(i-1)*length(T2)*length(phase_exc)+1:i*length(T2)*length(phase_exc)) = Dict;
        close(f)
    end
    
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

