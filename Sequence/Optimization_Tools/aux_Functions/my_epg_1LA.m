function [Fp,Fn,Z] = my_epg_1LA(exc_pulse,refoc_pulse,exc_phase,refoc_phase,T1,T2,dTE,ETL)

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
echoesFp  = zeros(ETL,length(exc_pulse)); % zero matrix ETL*Length exc_pulse
echoesFn  = zeros(ETL,length(exc_pulse)); % zero matrix ETL*Length exc_pulse
echoesZ   = zeros(ETL,length(exc_pulse)); % zero matrix ETL*Length exc_pulse
somaFp    = zeros(ETL,length(T2),length(exc_phase));
somaFn    = zeros(ETL,length(T2),length(exc_phase));
somaZ     = zeros(ETL,length(T2),length(exc_phase));


%% 2 - Generate epg dictionary using epg_cpmg.m function
for j=1:length(exc_phase)
    for i=1:length(T2)

        % -- Get value for dic. from CPMG condition
        for z=1:length(exc_pulse) % 1-->numero de TR (128)
            [s_Fp,s_Fn,s_Z] = epg_cpmg_1LA(exc_pulse(z),exc_phase(j),refoc_pulse(z,:),ETL,T1,T2(i),dTE,0,refoc_phase);
            echoesFp(:,z)   = s_Fp;
            echoesFn(:,z)   = s_Fn;
            echoesZ(:,z)    = s_Z;
        end        
        
        % -- Somatorio ao longo da fatia, soma das 128 coordenadas z        
        somaFp(:,i,j) = sum(echoesFp,2);  
        somaFn(:,i,j) = sum(echoesFn,2);  
        somaZ(:,i,j)  = sum(echoesZ,2);  
    end
end


%% 3 - Variable Dict orgnazied
D1_Fp = somaFp;
D1_Fn = somaFn;
D1_Z  = somaZ;

sz2_Fp = size(D1_Fp,2);
sz2_Fn = size(D1_Fn,2);
sz2_Z  = size(D1_Z,2);

sz3_Fp = size(D1_Fp,3);
sz3_Fn = size(D1_Fn,3);
sz3_Z  = size(D1_Z,3);

Fp = zeros(ETL,sz3_Fp*sz2_Fp);
Fn = zeros(ETL,sz3_Fn*sz2_Fn);
Z  = zeros(ETL,sz3_Z*sz2_Z);

for j = 1:sz2_Fp:size(Fp,2)
    Fp(:,j:j+sz2_Fp-1) = D1_Fp(:,:,(j+sz2_Fp-1)/sz2_Fp);
    Fn(:,j:j+sz2_Fn-1) = D1_Fn(:,:,(j+sz2_Fn-1)/sz2_Fn);
    Z(:,j:j+sz2_Z-1)   = D1_Z(:,:,(j+sz2_Z-1)/sz2_Z);
end

end
