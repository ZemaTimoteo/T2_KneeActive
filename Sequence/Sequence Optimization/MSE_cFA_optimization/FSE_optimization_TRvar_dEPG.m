
%% FSE optimization with CRLB analysis fitting EPG curves with dEPG
% Dictionaries previously generated
% TTFernandes - August 2024

% =============================
% Toolboxes - Add downloaded versions of 'Pulseq' and 'rf_tools' to the
%             toolboxes folder
% =============================

%% 0 - Set test
clear, clc

% Initial Sets 
testFSEopt  = 1;


%% 0 - Set matlab paths and toolboxes

% Active path
filePath    = matlab.desktop.editor.getActiveFilename;
justPath    = find(filePath == 'F');
file_path   = filePath([1:justPath(end)-1]); clear filePath
idx_tb_Path = find(file_path == 'M');
toolbox_Path = [file_path([1:idx_tb_Path(end)-1]) 'Toolboxes']; 
cd(file_path)


file_path_data      = [file_path 'Data'];
file_path_rf        = [file_path_data 'rf_pulses'];

filename_excell = 'parameters_&_bestResults.xlsx';

addpath(genpath(file_path))
addpath(genpath(toolbox_Path))


%% 1 - Tests

plotTest    = 'Fals';       % 'True' or 'Fals'
saveResults = 'True';

testSAR     = 'loadD';      % For B1_rms test = 'b1rms' OR Loading data 'loadD' 
testCRLB    = 'True';       % Test CRLB - 'True' OR Load CRLB - 'Fals'
testFBP     = 'Fals';       % Test Flip-back pulse - 'True' or 'Fals'
testTR      = 'Fals';

testFastSampling  = 'GRAPP';      % for SENSE - 'SENSE' OR for GRAPPA - 'GRAPP' OR for LORAKS - ' LORAK'
methodDic         = 'SLR_Prof';   % For CRLB - 'SLR_Prof' OR 'JUSTdict' OR '90DirSLR' OR '90SLRDir'

%% 2 - Parameters to Optimize

% ... 1.1 - Parameters for model optimization ...
B1 = 1;                             % B1 values for each we want to optimize FSE sequence
T1 = [1000 1000 1000];              % T1 max for knee and for T1 of GM @3T(ms)

% ... 1.2 - Parameters of sequence ...
nslices           = 25;             % number of slices - '30'
res               = 256;            % Resolution Nx = Ny
sliceThickn.exc   = 3e-3;           % Slice Thickness exc   (m) - '3e-3'
sliceThickn.refoc = 3e-3*3;         % Slice Thickness refoc (m) - '3e-3' * 3 

phase_exc     = pi/2;           % Phase Exciatation - in (rad) - along y
FA_exc_dic    = pi/2;           % Flip-angle Exciatation - in (rad) - along 

% ... 1.3 - Vectors of parameters ...
vector_dTE        = 8.35:0.02:12;       % Vector of Echo SpacingTime (ms) - parecido HSa - TODO acertar ao TE possível '8:2:30'
vector_ETL        = 6:1:20;             % Vector of # Echos - ver protocolo do HSa '6:2:30'
vector_flipAngle  = 130:1:180;          % Vector of flip Angles (ms) - ver artigo tammir - '90:5:180'
vector_SNR        = [40];               % variance - % variar SNR de teste  [1 30 100] OR [1 5:15:150];
vector_T2         = [ 8 45 ];           % vector of different T2 values to test - [8:8:50] OR '[8 45]'

% ... 1.4 - Parameters for generating EPG signals with slr_profile.py and my_epg.py ...
TR      = 10;                % units (ms) in Keerthivasan et al. 2019
T2_dic  = vector_T2;         % ms ( 40:1:55; OR 1:1:300; )
B1_dic  = 1;                 % B1 value ( 0.7:0.01:0.75; OR 0.6:0.01:1.4; )

% test - 22/11
if testTR == 'True'
    vector_TR = [1305 1939 2063 2189 2322 2460 2603 3779]; % in (ms)
end

% write in excell the parameters
% ----
cd(file_path_data)
A1 = {'Test','ini_TR (ms)','nslices','res (mm)','sliceThickn ex (mm)','sliceThickn ref (mm)','SNR','min_T2 (ms)','max_T2 (ms)','Test FBP'};
sheet = 1;
xlRange = 'A1';
xlswrite(filename_excell,A1,sheet,xlRange)

A2 = {testFSEopt,TR,nslices,res,sliceThickn.exc,sliceThickn.refoc,vector_SNR(1),min(vector_T2),max(vector_T2),testFBP};
sheet = 1;
xlRange = ['A',num2str(testFSEopt)];
xlswrite(filename_excell,A2,sheet,xlRange)
% ----

fprintf('\n\n 2 - Sucessfully finished -  Parameters to Optimize \n\n')


%% 3 - Constrains - (3min)
% Control the possible combinations of parameters
tic

% ... 3.1 - Parameters generation matrix
fline          = 1;
aline          = 1;
a              = 1;
failValues     = [];
acceptedValues = [];

 % »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»
% ... 3.2 - Constrains
 % »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»

    % 3.2.1 - SAR
SAR_max     = 100; % W/kg
RFpower_max = 50;  % W
maxB1_rms   = 5;   % uT - TODO tentar provocar os constrains

    % 3.2.2 - Acquisition Time
T_acq     = 12;             % Time in min - 10/12min
if testFastSampling == 'SENSE'
    af = 0.5;             % Acceleration factor
elseif testFastSampling == 'GRAPP'
    af = 2/3;             % Acceleration factor 
elseif testFastSampling == 'LORAK'
    accFactor = 4;                 % af of LORAKS
    af = 1/5 + (4/5)/accFactor;    % Acceleration sampling of k-space
end
maxTime   = T_acq;  % units(min)
maxTime_s = maxTime * 60;   % Time in (s)
    
    % 3.2.3 - Real resolution
real_res = res*af;

 % »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»
% ... 3.3 - Test for constrains | Test for B1+rms ...
 % »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»
 
% -> 3.3.1 - test for all vector ...
if testSAR == 'b1rms'
    clear b1_rms T_scan
    Models_failed = [];

    for jj=1:length(vector_dTE)

        for ii=1:length(vector_ETL)

            for gg=1:length(vector_flipAngle)

                % test - 22/11
                if testTR == 'True'
                    TR = vector_TR(gg); % in (ms)
                end

                if testFBP == 'Fals' % Test without Flip back pulse
                    % - Values for RFpower, SAR and T_scan ...
                    [b1_rms{jj,ii,gg},T_scan, TR_scan, Trec] = b1rms4seq_optFSE_TRvar(vector_dTE(jj), TR, ...
                        vector_ETL(ii), vector_flipAngle(gg), sliceThickn, nslices, real_res,'Fals',file_path_rf);
                else % Test with Flip back pulse
                    % - Values for RFpower, SAR and T_scan ...
                    [b1_rms{jj,ii,gg},T_scan, TR_scan, Trec] = b1rms4seq_optFSE_TRvar_FBP(vector_dTE(jj), TR, ...
                        vector_ETL(ii), vector_flipAngle(gg), sliceThickn, nslices, real_res,'Fals');
                end

                % - check B1+rms values and max acq Time ...
                if isnan(b1_rms{jj,ii,gg}) %Models_failed due to constrain b1_rms
                    Models_failed(fline,:) = [vector_dTE(jj) vector_ETL(ii) vector_flipAngle(gg) ...
                        res T_scan/60 NaN Trec TR_scan];
                    fline = fline+1;
                    continue
                elseif b1_rms{jj,ii,gg} > maxB1_rms || T_scan > maxTime_s %Models_failed due to constrain b1_rms
                    Models_failed(fline,:) = [vector_dTE(jj) vector_ETL(ii) vector_flipAngle(gg) ...
                        res T_scan/60 b1_rms{jj,ii,gg} Trec TR_scan];
                    fline = fline+1;
                else %Models_accepted: [TE(ms) #ETL FA Nx/Ny TimeScan(min) b1+rms(uT) Trec(ms) TR_scan(s)]
                    Models_accepted(aline,:) = [vector_dTE(jj) vector_ETL(ii) vector_flipAngle(gg) ...
                        res T_scan/60 b1_rms{jj,ii,gg} Trec TR_scan];
                    aline = aline+1;
                end
                a=a+1;
            end

        end
        numModels = aline + fline - 2; % Number of all models
        fprintf(['      Constrain check complete ',num2str(jj),' / ',num2str(length(vector_dTE)),'\n'])
    end

     % »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»
    % ... 3.4 - Save data ...
     % »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»

    if saveResults == 'True'
        cd(file_path_data)

        if testFBP == 'Fals'
            save(['Test',num2str(testFSEopt),'_ConstrainModels_maxB1rms',num2str(maxB1_rms),'_maxTime',num2str(maxTime_s),...
                '_mindTE',num2str(min(vector_dTE)),'_maxdTE',num2str(max(vector_dTE)), ...
                '_minETL',num2str(min(vector_ETL)),'_maxETL',num2str(max(vector_ETL)), ...
                '_minFlipA',num2str(min(vector_flipAngle)),'_maxFlipA',num2str(max(vector_flipAngle)),...
                '.mat'],'Models_accepted','Models_failed','maxB1_rms','maxTime_s')
        else
            save(['Test',num2str(testFSEopt),'FPB_ConstrainModels_maxB1rms',num2str(maxB1_rms),'_maxTime',num2str(maxTime_s),...
                '_mindTE',num2str(min(vector_dTE)),'_maxdTE',num2str(max(vector_dTE)), ...
                '_minETL',num2str(min(vector_ETL)),'_maxETL',num2str(max(vector_ETL)),...
                '_minFlipA',num2str(min(vector_flipAngle)),'_maxFlipA',num2str(max(vector_flipAngle)),...
                '.mat'],'Models_accepted','Models_failed','maxB1_rms','maxTime_s')
        end
        cd(file_path)
    end



elseif testSAR == 'loadD'
    
    cd(file_path_data)
    if testFBP == 'Fals'
        load(['Test',num2str(testFSEopt),'_ConstrainModels_maxB1rms',num2str(maxB1_rms),'_maxTime',num2str(maxTime_s),...
            '_mindTE',num2str(min(vector_dTE)),'_maxdTE',num2str(max(vector_dTE)), ...
            '_minETL',num2str(min(vector_ETL)),'_maxETL',num2str(max(vector_ETL)),...
            '_minFlipA',num2str(min(vector_flipAngle)),'_maxFlipA',num2str(max(vector_flipAngle)),...
            '.mat'])
    else
        load(['Test',num2str(testFSEopt),'FPB_ConstrainModels_maxB1rms',num2str(maxB1_rms),'_maxTime',num2str(maxTime_s),...
            '_mindTE',num2str(min(vector_dTE)),'_maxdTE',num2str(max(vector_dTE)), ...
            '_minETL',num2str(min(vector_ETL)),'_maxETL',num2str(max(vector_ETL)),...
            '_minFlipA',num2str(min(vector_flipAngle)),'_maxFlipA',num2str(max(vector_flipAngle)),...
            '.mat'])
    end
    cd(file_path)   
end   
toc

%Models_accepted: [TE(ms) #ETL angle TimeScan(min) b1+rms(uT) Trec(ms)]
fprintf(['      ... Models Accepted | Models Rejected ',num2str(size(Models_accepted,1)),' | ',num2str(size(Models_failed,1)),'\n'])
fprintf('\n\n 3 - Sucessfully finished -  Constrains \n\n')


%% 4 - CRLB (8min06s)

if testCRLB == 'True'
    tic
           
    % ... 4.1 - Parameters ...
    maxETL    = max(Models_accepted(:,2));
    sigma     = zeros(1,size(vector_SNR,2));
    dEPG_dT2  = zeros(size(vector_ETL,2),size(vector_SNR,2));
    vardT2    = zeros(size(Models_accepted,1),size(vector_SNR,2),size(vector_T2,2));
    factorFBP = zeros(size(Models_accepted,1),size(vector_SNR,2),size(vector_T2,2));      
    EPG       = zeros(1,size(vector_SNR,2));
    
    % ... 4.2 - Cycle ... 
    %Models_accepted: [ 1-TE(ms)  2-#ETL 3-angle 4-res(mm) 5-TimeScan(min) 6-b1+rms 7-TRec(ms) 8-TR(s) ]        
    for gg=1:size(vector_T2,2) % cycle over T2

        for jj =1:size(vector_SNR,2) % Itteration by SNR values

            for ii=1:size(Models_accepted,1) % cycle over TE
                              
                % ... 4.2.1 - Parameters ...
                FA_refoc_dic  = ones(1,Models_accepted(ii,2))*Models_accepted(ii,3);   % Flip-angle Refocosuing - in (degrees) - along y
                TRacq         = Models_accepted(ii,8);                                 % TRacq - in (s)
                sigma(jj)     = (1/vector_SNR(jj));                                    % Noise variance for CRLB 
                plotTest      = 'Fals';
                
                % ... 4.2.2 - Get CRLB num and uncertainty  ...
                %function [vardT2, ds_dT2, FF_tool, factorFPB_epg] = CRLB_epg_optFSE_TRvar_vF(rf_exc, B1, T1, T2, dTE, ETL, TRacq, Trec, flipA, sigma, plotTest,testFPB, dir_data, methodDic)                
                [vardT2(ii,jj,gg),ds_dT2{ii,jj,gg},ff{ii,jj,gg},factorFBP(ii,jj,gg)] = ...
                            CRLB_epg_optFSE_TRvar_vF(...
                                            FA_exc_dic, B1, T1(gg),  vector_T2(gg),...
                                            Models_accepted(ii,1), Models_accepted(ii,2), ...
                                            TRacq, Models_accepted(ii,7), FA_refoc_dic, sigma(jj), ...
                                            plotTest, testFBP, file_path_data, methodDic...
                                            );
                
                if plotTest  == 'True'
                    figure,plot(abs(ff{ii,jj,gg}),'*--'),hold on,
                    title('Signal - EPG')
                    legend('Original Signal')

                    figure()
                    plot(abs(ds_dT2{ii,jj,gg})/max(abs(ds_dT2{ii,jj,gg})),'ro--'),
                    hold on
                    title('Analitical & Numberical for T2 derivatives - dEPG')
                    legend('T2 derivative dEPG')

                end                

                fprintf(['      Test ',num2str(ii),' / ',num2str(size(Models_accepted,1)),'\n'])
            end
            
            fprintf(['  --  SNR Test ',num2str(jj),' / ',num2str(size(vector_SNR,2)),'\n'])
        end
        fprintf(['\n\n      Successfull CRLB model obtained for different T2 values ',num2str(gg),...
            ' / ',num2str(size(vector_T2,2)),'\n'])
    end
    toc
    
    % save CRLB data
    cd(file_path_data)
    if testFBP == 'Fals'
        save(['CRLBvalues_test',num2str(testFSEopt),...
            '_TRvar_ModelsAccepted',num2str(size(Models_accepted,1)),...
            '_minSNR',num2str(min(vector_SNR)),'_maxSNR',num2str(max(vector_SNR)), ...
            '_minT2-',num2str(min(vector_T2)),'_maxT2-',num2str(max(vector_T2)),...
            '.mat'],'vardT2','dEPG_dT2','EPG','factorFBP')
    else
        save(['CRLBvalues_FPB_test',num2str(testFSEopt),...
            '_TRvar_ModelsAccepted',num2str(size(Models_accepted,1)),...
            '_minSNR',num2str(min(vector_SNR)),'_maxSNR',num2str(max(vector_SNR)), ...
            '_minT2-',num2str(min(vector_T2)),'_maxT2-',num2str(max(vector_T2)),...
            '.mat'],'vardT2','dEPG_dT2','EPG','factorFBP')
    end
    cd(file_path)
    
elseif testCRLB == 'Fals'
    % load CRLB data    
    cd(file_path_data)
    if testFBP == 'Fals'      
        load(['CRLBvalues_test',num2str(testFSEopt),...
            '_TRvar_ModelsAccepted',num2str(size(Models_accepted,1)),...
            '_minSNR',num2str(min(vector_SNR)),'_maxSNR',num2str(max(vector_SNR)), ...
            '_minT2-',num2str(min(vector_T2)),'_maxT2-',num2str(max(vector_T2)),...
            '.mat'])
    else
        load(['CRLBvalues_FPB_test',num2str(testFSEopt),...
            '_TRvar_ModelsAccepted',num2str(size(Models_accepted,1)),...
            '_minSNR',num2str(min(vector_SNR)),'_maxSNR',num2str(max(vector_SNR)), ...
            '_minT2-',num2str(min(vector_T2)),'_maxT2-',num2str(max(vector_T2)),...
            '.mat'])
    end
    cd(file_path)
end

if testFBP == 'True'
    aa = reshape(factorFBP,size(factorFBP,1)*size(factorFBP,3),1);
    figure();plot(aa) 
end



if plotTest == 'True'
    figure()
    subplot(121)    
    plot(squeeze(vardT2(:,1,1)))
    hold on
    title('T2 = 45 EPG')
    subplot(122)
    plot(squeeze(vardT2(:,1,2)))
    hold on
    title('T2 = 200 EPG')    
end


fprintf('\n\n 4 - Sucessfully finished - CRLB uncertainty\n\n')





%% 5 - Get best model fit
clear aux_results_dTE aux_results_ETL aux_results_flipA a b c uni_dTE_x uni_dTE_y ...
    uni_ETL_x uni_ETL_y uni_flipA_x uni_flipA_y results_dTE results_ETL results_flipA ...
    Results
    
% ... 5.1 - Parameters ...
SNR_val   = vector_SNR(1);    % Possible SNR values = 1 30 100
T2_val    = vector_T2(1);     % units (ms) - Possible T2 values = 8    16    24    32    40    48

testdTE   = vector_dTE(4);           % units (ms)
testETL   = vector_ETL(3);           % echoes
testFlipA = vector_flipAngle(3);    % units (º) - 115

SNRindx    = find(vector_SNR == SNR_val);
aux_T2indx = find(vector_T2 == T2_val);
T2indx     = (size(Models_accepted,2)+aux_T2indx);
% T2thrs     = 10*T2_val;

% ... 5.1 - Get  Results ...
%Models_accepted: [  TE(ms)  #ETL   angle  res  TimeScan(min)   b1+rms(uT)  Trec(ms) TR]
%Models_accepted: [       1     2       3    4              5            6         7  8]
% % Results         = Models_accepted; 
Results_vardT2  = Models_accepted; 

% % for gg=1:size(vector_T2,2)
% %     Results  = [Results uCRLB(:,SNRindx,gg)];
% % end

% time balanced CRLB by scan time!
for gg=1:size(vector_T2,2)
    for ii=1:size(Models_accepted,1)
% %         new_uCRLB(ii,SNRindx,gg)  = uCRLB(ii,SNRindx,gg)*sqrt(Models_accepted(ii,5)*60); 
        new_vardT2(ii,SNRindx,gg) = vardT2(ii,SNRindx,gg)*sqrt(Models_accepted(ii,5)*60); % ??
    end
end

%Results_vardT2: [  TE(ms)  #ETL   angle  res  TimeScan(min)   b1+rms(uT)  Trec(ms)  TR(s)  var results for T2(1)   T2(2)   ]
%                [   '1'     '2'    '3'   '4'       '5'           '6'        '7'      '8'                    '9'    '10'    ]
for gg=1:size(vector_T2,2)
    Results_vardT2  = [Results_vardT2 new_vardT2(:,SNRindx,gg)];
end

% aux for plot
% % intervPlot        = [min(min(log(Results(:,T2indx:end)))) max(max(log(Results(:,T2indx:end))))]; % for caxis
% % intervPlot_vardT2 = [0 max(max(log(Results_vardT2(:,T2indx:end))))]; % for caxis
intervPlot_vardT2        = [min(min(log(Results_vardT2(:,T2indx:end)))) max(max(log(Results_vardT2(:,T2indx:end))))]; % for caxis

% ... 5.2 - Get matrix for plot ...
a=1;b=1;c=1;
for ii=1:size(Models_accepted,1)    
    % Test for fixed Echo Time (dTE)
    if Results_vardT2(ii,1) == testdTE  % Test for fixed dTE
        aux_results_dTE(a,:) = Results_vardT2(ii,[2 3 T2indx]);
        a = a+1;
    end
    % Test for fixed Echo train Length (ETL)
    if Results_vardT2(ii,2) == testETL
        aux_results_ETL(b,:) = Results_vardT2(ii,[1 3 T2indx]);
        b = b+1;
    end
    % Test for fixed flip Angle
    if Results_vardT2(ii,3) == testFlipA
        aux_results_flipA(c,:) = Results_vardT2(ii,[1 2 T2indx]);
        c = c+1;
    end    
end

% ... 5.3 - Get matrix for plot ...
uni_dTE_x   = unique(aux_results_dTE(:,1));
uni_dTE_y   = unique(aux_results_dTE(:,2));
uni_ETL_x   = unique(aux_results_ETL(:,1));
uni_ETL_y   = unique(aux_results_ETL(:,2));
uni_flipA_x = unique(aux_results_flipA(:,1));
uni_flipA_y = unique(aux_results_flipA(:,2));

results_dTE   = NaN(length(uni_dTE_x),length(uni_dTE_y));
results_ETL   = NaN(length(uni_ETL_x),length(uni_ETL_y));
results_flipA = NaN(length(uni_flipA_x),length(uni_flipA_y));

for i=1:size(aux_results_dTE,1)
    for g=1:length(uni_dTE_x)
        for j=1:length(uni_dTE_y)
            if uni_dTE_x(g)==aux_results_dTE(i,1) && uni_dTE_y(j)==aux_results_dTE(i,2)
                results_dTE(g,j) = aux_results_dTE(i,3);
            end
        end
    end
end
for i=1:length(aux_results_ETL)
    for g=1:length(uni_ETL_x)
        for j=1:length(uni_ETL_y)
            if uni_ETL_x(g)==aux_results_ETL(i,1) && uni_ETL_y(j)==aux_results_ETL(i,2)
                results_ETL(g,j) = aux_results_ETL(i,3);
            end
        end
    end
end
for i=1:size(aux_results_flipA,1)
    for g=1:length(uni_flipA_x)
        for j=1:length(uni_flipA_y)
            if uni_flipA_x(g)==aux_results_flipA(i,1) && uni_flipA_y(j)==aux_results_flipA(i,2)
                results_flipA(g,j) = aux_results_flipA(i,3);
            end
        end
    end
end
clear  aux_results_dTE aux_results_ETL aux_results_flipA


% ... 5.4 - Plots ...
figure()
% plot for fixed dTE
subplot(221)
pcolor(log(results_dTE))
hold on
grid on
title(['log Var/T2 | T2 = ',num2str(vector_T2(aux_T2indx)),' ms | fixed dTE = ',num2str(testdTE),' ms'],'fontsize',20)
xlabel(['flipAngle (º)'],'fontsize',20)
ylabel(['# Echoes'],'fontsize',20)
set(gca,'FontSize',15)
xticks(1:1:size(uni_dTE_y,1)*2)
yticks(1:1:size(uni_dTE_x,1)*2)
xticklabels(uni_dTE_y)
yticklabels(uni_dTE_x)
caxis(intervPlot_vardT2)

% plot for fixeed ETL
subplot(222)
pcolor(log(results_ETL))
hold on
grid on
title(['log Var/T2 | T2 = ',num2str(vector_T2(aux_T2indx)),' ms | fixed ETL = ',num2str(testETL)],'fontsize',20)
xlabel(['flipAngle (º)'],'fontsize',20)
ylabel(['dTE (ms)'],'fontsize',20)
set(gca,'FontSize',15)
xticks(1:1:size(uni_ETL_y,1)*2)
yticks(1:1:size(uni_ETL_x,1)*2)
xticklabels(uni_ETL_y)
yticklabels(uni_ETL_x)
caxis(intervPlot_vardT2)

% plot for fixeed flip Angle
subplot(223)
pcolor(log(results_flipA))
hold on
grid on
title(['log Var/T2 | T2 = ',num2str(vector_T2(aux_T2indx)),'ms | fixed fA = ',num2str(testFlipA),'º'],'fontsize',20)
xlabel(['# Echoes'],'fontsize',20)
ylabel(['dTE (ms)'],'fontsize',20)
set(gca,'FontSize',15)
xticks(1:1:size(uni_flipA_y,1)*2)
yticks(1:1:size(uni_flipA_x,1)*2)
xticklabels(uni_flipA_y)
yticklabels(uni_flipA_x)
caxis(intervPlot_vardT2)


%% 6 - Values of maximized T2 variance

%SNR = 30 for different T2 values
cd(file_path_data)
fprintf(['\n\n            ... Optimized Results with CRLB variance/T2 ... \n\n'])

for gg=1:size(vector_T2,2)
    % write in excell the 
    
    B = {'Test','log(T2 max variance)','T2 (ms)','dTE (ms)','ETL','SNR','flipAngle (º)','Trec (s)','TR_acq (s)'};
    sheet = 1+gg;
    xlRange_title = 'A1';
    xlswrite(filename_excell,B,sheet,xlRange_title)
    
    fprintf(['-------------------\n\n'])
    
    [maxValue_vardT2(gg), indxResult_vardT2] = (max( (Results_vardT2(:,T2indx-aux_T2indx+gg)) ) ) ;
    fprintf(['   T2 max variance ' ,num2str((maxValue_vardT2(gg))),' for: ', ...
        'T2 = ', num2str(vector_T2(gg)),...        
        ' ms, dTE = ', num2str(Results_vardT2(indxResult_vardT2,1)),...
        ' ms, ETL = ', num2str(Results_vardT2(indxResult_vardT2,2)),...
        ' SNR = ',num2str(SNR_val),...
        ', flipAngle = ', num2str(Results_vardT2(indxResult_vardT2,3)),'º',...
        ', Trec = ', num2str(Results_vardT2(indxResult_vardT2,7)*1e-3),'s',...
        ', TR = ',num2str(Results_vardT2(indxResult_vardT2,8)),'s\n\n',...
        ', Tscan = ',num2str(Results_vardT2(indxResult_vardT2,8)*real_res/60)])
    
    results.dTE_vector(gg)       =  Results_vardT2(indxResult_vardT2,1);
    results.ETL_vector(gg)       =  Results_vardT2(indxResult_vardT2,2);
    results.SNR_val_vector(gg)   =  SNR_val;
    results.flipAngle_vector(gg) =  Results_vardT2(indxResult_vardT2,3);
    results.TR_acq_vector(gg)    =  Results_vardT2(indxResult_vardT2,7);
    results.TR_acq_vector(gg)    =  Results_vardT2(indxResult_vardT2,8);
    results.vardT2               =  Results_vardT2(indxResult_vardT2,8);

    fprintf(['-------------------\n'])

    % write in excell the parameters
    C = {testFSEopt,log(maxValue_vardT2(gg)),vector_T2(gg),Results_vardT2(indxResult_vardT2,1), ...
            Results_vardT2(indxResult_vardT2,2),SNR_val,...
            Results_vardT2(indxResult_vardT2,3),...
            Results_vardT2(indxResult_vardT2,9)};
    xlRange_T2 = ['A',num2str(testFSEopt)];
    xlswrite(filename_excell,C,sheet,xlRange_T2)    
end


%% 7 - Figures (ISMRM23 - Plots)
figure()

for aa= 1:size(vector_T2,2)
    clear aux_results_dTE aux_results_ETL aux_results_flipA a b c uni_dTE_x uni_dTE_y ...
        uni_ETL_x uni_ETL_y uni_flipA_x uni_flipA_y results_dTE results_ETL results_flipA ...
        Results

    % ... 7.1 - Parameters ...
    SNR_val   = results.SNR_val_vector(aa);      % Possible SNR values = 1 30 100
    T2_val    = vector_T2(aa);                   % units (ms) - Possible T2 values = 8    16    24    32    40    48
    testdTE   = results.dTE_vector(aa);          % units (ms)
    testETL   = results.ETL_vector(aa);          % echoes
    testFlipA = results.flipAngle_vector(aa);    % units (º) - 115

    SNRindx    = find(vector_SNR == SNR_val);
    aux_T2indx = find(vector_T2 == T2_val);
    T2indx     = (size(Models_accepted,2)+aux_T2indx);
    % T2thrs     = 10*T2_val;

    % ... 7.1 - Get  Results ...
    %Models_accepted: [  TE(ms)  #ETL   angle  res  TimeScan(min)   b1+rms(uT)  Trec(ms) TR]
    % % Results         = Models_accepted; 
    Results_vardT2  = Models_accepted; 

    % % for gg=1:size(vector_T2,2)
    % %     Results  = [Results uCRLB(:,SNRindx,gg)];
    % % end

    % time balanced CRLB by scan time!
    for gg=1:size(vector_T2,2)
        for ii=1:size(Models_accepted,1)
    % %         new_uCRLB(ii,SNRindx,gg)  = uCRLB(ii,SNRindx,gg)*sqrt(Models_accepted(ii,5)*60); 
            new_vardT2(ii,SNRindx,gg) = vardT2(ii,SNRindx,gg)*sqrt(Models_accepted(ii,5)*60); % ??
        end
    end

    %Results_vardT2: [  TE(ms)  #ETL   angle  res  TimeScan(min)   b1+rms(uT)  Trec(ms)  TR(s)  results for T2(1)   T2(2)   ]
    %                [   '1'     '2'    '3'   '4'       '5'           '6'        '7'      '8'                 '9'    '10'    ]
    for gg=1:size(vector_T2,2)
    % %     Results         = [Results new_uCRLB(:,SNRindx,gg)];
        Results_vardT2  = [Results_vardT2 new_vardT2(:,SNRindx,gg)];
    end

    % aux for plot
    % % intervPlot        = [min(min(log(Results(:,T2indx:end)))) max(max(log(Results(:,T2indx:end))))]; % for caxis

    % ... 7.2 - Get matrix for plot ...
    a=1;b=1;c=1;
    for ii=1:size(Models_accepted,1)    
        % Test for fixed Echo Time (dTE)
        if Results_vardT2(ii,1) == testdTE  % Test for fixed dTE
            aux_results_dTE(a,:) = Results_vardT2(ii,[2 3 T2indx]);
            a = a+1;
        end
        % Test for fixed Echo train Length (ETL)
        if Results_vardT2(ii,2) == testETL
            aux_results_ETL(b,:) = Results_vardT2(ii,[1 3 T2indx]);
            b = b+1;
        end
        % Test for fixed flip Angle
        if Results_vardT2(ii,3) == testFlipA
            aux_results_flipA(c,:) = Results_vardT2(ii,[1 2 T2indx]);
            c = c+1;
        end    
    end

    % ... 7.3 - Get matrix for plot ...
    uni_dTE_x   = unique(aux_results_dTE(:,1));
    uni_dTE_y   = unique(aux_results_dTE(:,2));
    uni_ETL_x   = unique(aux_results_ETL(:,1));
    uni_ETL_y   = unique(aux_results_ETL(:,2));
    uni_flipA_x = unique(aux_results_flipA(:,1));
    uni_flipA_y = unique(aux_results_flipA(:,2));

    results_dTE   = NaN(length(uni_dTE_x),length(uni_dTE_y));
    results_ETL   = NaN(length(uni_ETL_x),length(uni_ETL_y));
    results_flipA = NaN(length(uni_flipA_x),length(uni_flipA_y));

    for i=1:size(aux_results_dTE,1)
        for g=1:length(uni_dTE_x)
            for j=1:length(uni_dTE_y)
                if uni_dTE_x(g)==aux_results_dTE(i,1) && uni_dTE_y(j)==aux_results_dTE(i,2)
                    results_dTE(g,j) = aux_results_dTE(i,3);
                end
            end
        end
    end
    for i=1:length(aux_results_ETL)
        for g=1:length(uni_ETL_x)
            for j=1:length(uni_ETL_y)
                if uni_ETL_x(g)==aux_results_ETL(i,1) && uni_ETL_y(j)==aux_results_ETL(i,2)
                    results_ETL(g,j) = aux_results_ETL(i,3);
                end
            end
        end
    end
    for i=1:size(aux_results_flipA,1)
        for g=1:length(uni_flipA_x)
            for j=1:length(uni_flipA_y)
                if uni_flipA_x(g)==aux_results_flipA(i,1) && uni_flipA_y(j)==aux_results_flipA(i,2)
                    results_flipA(g,j) = aux_results_flipA(i,3);
                end
            end
        end
    end
    clear  aux_results_dTE aux_results_ETL aux_results_flipA


    % ... 7.4 - Plots ...
    % plot for fixed dTE
    subplot(size(vector_T2,2),1,aa)
    pcolor(log(results_dTE))
%     pcolor((results_dTE))
    hold on
    grid on
    title(['log Var/T2 | T2 = ',num2str(vector_T2(aux_T2indx)),' ms | fixed dTE = ',num2str(testdTE),' ms'],'fontsize',20)
    xlabel(['flipAngle (º)'],'fontsize',20)
    ylabel(['# Echoes'],'fontsize',20)
    set(gca,'FontSize',15)
    xticks(1:1:size(uni_dTE_y,1)*2)
    yticks(1:1:size(uni_dTE_x,1)*2)
    xticklabels(uni_dTE_y)
    yticklabels(uni_dTE_x)
    hold on

end
% intervPlotISMRM23 = [0 max(max(log(results_dTE(:,:))))]; % for caxis
intervPlotISMRM23 = [min(min((results_dTE(:,:)))) max(max((results_dTE(:,:))))]; % for caxis
caxis(intervPlotISMRM23)

