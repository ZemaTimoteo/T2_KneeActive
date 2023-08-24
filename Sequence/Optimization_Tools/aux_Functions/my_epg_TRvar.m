function Dict = my_epg_TRvar(exc_pulse,refoc_pulse,exc_phase,refoc_phase,T1,T2,dTE,ETL,Trec)

%   Input:
%       exc_pulse   -> Magnitude/Profile of excitatory pulse
%       refoc_pulse -> Magnitude/Profile of Refocus Pulse
%       exc_phase   -> Phase of Excitatory Pulse
%       refoc_phase -> Phase of Refocus Pulse
%       T1          -> T1 value
%       T2          -> T2 value
%       dTE         -> Echo Spacing in ms
%       ETL         -> Echo Train Lenght
%
%   Output:
%       Dict        -> Dictionary simulation values


%% 1 - Set parameters
exc_pulse = exc_pulse(:);
echoes    = zeros(ETL,length(exc_pulse)); % zero matrix ETL*Length exc_pulse
soma      = zeros(ETL,length(T2),length(exc_phase));
auxVariab = 1 - exp(-Trec/T1);

%% 2 - Generate epg dictionary using epg_cpmg.m function
for j=1:length(exc_phase)
    for i=1:length(T2)
% %         disp('T2: '),i

        % -- Get value for dic. from CPMG condition
        for z=1:length(exc_pulse) % 1-->pontos dos pulsos (128)
            s = epg_cpmg(exc_pulse(z),exc_phase(j),refoc_pulse(z,:),ETL,T1,T2(i),dTE,0,refoc_phase);
% %             s1= epg_cpmg(pi          ,30          ,1000            ,200,10);	% IDEAL, CPMG
            echoes(:,z) = s;
        end
        
        % ------------ Slice thickness --------------------------
        
        % ... testes ...
        % xvec_shortest = 0.8656;  (slr_profile.m)
        % figure, plot(linspace(-xvec_shortest,xvec_shortest,128),abs(echoes(1,:)))        
        % figure, plot(abs(echoes(1,:))), hold on, plot(abs(echoes(2,:)),'r')
        % hold on, plot(abs(echoes(5,:)),'k')
        
        % -- Somatorio ao longo da fatia, soma das 128 coordenadas z        
        soma(:,i,j) = sum(echoes,2);
        
% %         figure; plot(abs(soma(:,:)))
% %         legend('Signal of EPG')

    end
end


%% 3 - Variable Dict orgnazied
D1 = soma;
sz2 = size(D1,2);
sz3 = size(D1,3);
Dict = zeros(ETL,sz3*sz2);

for j = 1:sz2:size(Dict,2)
    Dict(:,j:j+sz2-1) = D1(:,:,(j+sz2-1)/sz2);
end

end
