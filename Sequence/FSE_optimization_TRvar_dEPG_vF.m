%% FSE optimization with CRLB analysis fitting EPG curves with dEPG
% Dictionaries previously generated
% TTFernandes - May2022

% =============================
% Toolboxes:
%      - sar4seq - in: https://github.com/imr-framework/sar4seq/tree/master/sar4seq
%      - T2estimation-master 
% =============================
%% 0 - Set test
clear, clc

% Initial Sets 
testFSEopt  = 21;


%% 0 - Set matlab paths and toolboxes

file_path    = 'Github\Sequence\MSE_files';
%toolboxes
addpath(genpath('Github\Toolboxes\sar4seq-master'))
addpath(genpath('Github\Toolboxes\pulseq-master'));
addpath(genpath('Github\Sequence\Optimization_Tools\Optimization_Tools'));
addpath(genpath('Github\Toolboxes\T2estimation-master\src'))

file_path_dic       = [file_path '\Dictionaries'];
file_path_data      = [file_path '\Data'];
file_path_dic_save  = [file_path '\Dictionaries\FAVar_test' num2str(testFSEopt)];

filename_excell = 'parameters_&_bestResults.xlsx';

addpath(genpath(file_path))


%% 1 - Tests

plotTest    = 'True';       % 'True' or 'Fals'
saveResults = 'Fals';

testSAR     = 'b1rms';      % For B1_rms test = 'b1rms' OR For measuring SAR = 'SAR_t' OR Loading data 'loadD' 
testCRLB    = 'True';       % Test CRLB - 'True' OR Load CRLB - 'Fals'
testFBP     = 'Fals';       % Test Flip-back pulse - 'True' or 'Fals'

testSENSE   = 'True';
testGRAPPA  = 'Fals';

testTR      = 'Fals';

methodDic   = 'JUSTdict';   % For CRLB - 'SLR_Prof' OR 'JUSTdict' OR '90DirSLR' OR '90SLRDir'

%% 2 - Parameters to Optimize

% ... 1.1 - Parameters for model optimization ...
B1            = 1;     % B1 values for each we want to optimize FSE sequence
T1            = 1;  % T1 max for knee (s)
% % T1            = 0.598;  % T1 max for knee (s)
T1maxKnee     = T1;  % T1 max for knee (s)
Sample_weight = 75;             % Weight of subject (kg)

% ... 1.2 - Parameters of sequence ...
nslices       = 25;             % number of slices - '30'
res           = 256;            % Resolution Nx = Ny
sliceThickn   = 2.5e-3;         % Slice Thickness (m) - '2.6e-3'

phase_exc     = pi/2;           % Phase Exciatation - in (rad) - along y
FA_exc_dic    = pi/2;           % Flip-angle Exciatation - in (rad) - along 

% ... 1.3 - Vectors of parameters ...
vector_dTE        = 8:1:25;      % Vector of Echo SpacingTime (ms) - parecido HSa - TODO acertar ao TE possível '8:2:30'
vector_ETL        = 6:2:40;      % Vector of # Echos - ver protocolo do HSa '6:2:30'
vector_flipAngle  = 90:5:180;         % Vector of flip Angles (ms) - ver artigo tammir - '90:5:180'
vector_SNR        = [40];        % variance - % variar SNR de teste  [1 30 100] OR [1 5:15:150];
vector_T2         = [8 45 200];    % vector of different T2 values to test - [8:8:50] OR '[8 45]'

% ... 1.4 - Parameters for generating EPG signals with slr_profile.py and my_epg.py ...
% % TR      = 5*max(T1maxKnee); % units (ms) in Keerthivasan et al. 2019
% % TR      = 5*100; % units (ms) in Keerthivasan et al. 2019
TR      = 2;                % units (ms) in Keerthivasan et al. 2019
% TR      = 7454;             % Repetition time in (ms)T1_dic  = T1maxKnee;        % ms
T2_dic  = vector_T2;        % ms ( 40:1:55; OR 1:1:300; )
B1_dic  = 1;                % B1 value ( 0.7:0.01:0.75; OR 0.6:0.01:1.4; )

% test - 22/11
if testTR == 'True'
    vector_TR = [1305 1939 2063 2189 2322 2460 2603 3779]; % in (ms)
end

% write in excell the parameters
% ----
cd(file_path_data)
A1 = {'Test','ini_TR (ms)','nslices','res (mm)','sliceThickn (mm)','SNR','min_T2 (ms)','max_T2 (ms)','Test FBP'};
sheet = 1;
xlRange = 'A1';
xlswrite(filename_excell,A1,sheet,xlRange)

A2 = {testFSEopt,TR,nslices,res,sliceThickn,vector_SNR(1),min(vector_T2),max(vector_T2),testFBP};
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

% ... 3.2 - Constrains
    % 3.2.1 - SAR
SAR_max     = 100; % W/kg
RFpower_max = 50;  % W
maxB1_rms   = 2.8;   % uT - TODO tentar provocar os constrains

% 3.2.1 - Acquisition Time
if testSENSE == 'True'
    T_acq     = 12;            % Time in min - 10/12min
    SENSE     = 2;             % Acceleration factor
    maxTime   = T_acq*SENSE;   % units(min)
    maxTime_s = maxTime * 60;  % Time in (s)
elseif testGRAPPA == 'True'
    T_acq     = 12;            % Time in min - 10/12min
    af        = 2;             % Acceleration factor
    maxTime   = T_acq*2/3;   % units(min)
    maxTime_s = maxTime * 60;  % Time in (s)    
end

% ... 3.3 - Test for constrains ...
% % vector_flipAngle = ones(1,ETL)*flipAngle;

 % »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»
    % 3.3.1 - Test for SAR
 % »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»
if testSAR == 'SAR_t'
    
        % -> 3.3.1 a) - test for all vector ...
    [RFpower,SAR]    = SAR4seq_optFSE(Sample_weight,dTE,TR,ETL,vector_flipAngle,sliceThickn,nslices,res);

        % -> 3.3.1 b) - test for all vector ...
    for jj=1:length(vector_dTE)
        for ii=1:length(vector_ETL)
            f = waitbar((a)/(length(vector_dTE)*length(vector_ETL)),'SAR & Time constrains...');

            vector_flipAngle = ones(1,vector_ETL(ii))*flipAngle;

            % - Values for RFpower, SAR and T_scan ...
            [RFpower{jj,ii},SAR{jj,ii},T_scan] = SAR4seq_optFSE(Sample_weight,vector_dTE(jj), ...
                vector_ETL(ii),vector_flipAngle,sliceThickn, ...
                nslices,res);

            % - check SAR values and max acq Time ...
            close(f)
            if isnan(RFpower{jj,ii}.wbg_tavg)
                failValues(fline,:) = [jj ii];
                fline = fline+1;
                continue
            elseif RFpower{jj,ii}.wbg_tavg > RFpower_max || SAR{jj,ii}.SAR_10s_wbg > SAR_max || T_scan > maxTime_s
                failValues(fline,:) = [jj ii];
                fline = fline+1;
            else
                acceptedModels(aline,:) = [vector_dTE(jj) vector_ETL(ii) T_scan/60]; %[TE #ETL TimeScan]
                aline = aline+1;
            end
            a=a+1;
        end
    end

    
 % »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»
    % 3.3.2 - Test for B1+rms
 % »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»
elseif testSAR == 'b1rms'
    % -> 3.3.2 a) - test for all vector ...
% %     [b1_rms,T_scan]    = b1rms4seq_optFSE(Sample_weight,dTE,TR,ETL,flipAngle,sliceThickn,nslices,res,'Fals');
    
    % -> 3.3.2 b) - test for all vector ...
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
                    [b1_rms{jj,ii,gg},T_scan, TR_scan, Trec] = b1rms4seq_optFSE_TRvar(Sample_weight,vector_dTE(jj), TR, ...
                        vector_ETL(ii), vector_flipAngle(gg), sliceThickn, nslices, res,'Fals');
                else % Test with Flip back pulse
                    % - Values for RFpower, SAR and T_scan ...
                    [b1_rms{jj,ii,gg},T_scan, TR_scan, Trec] = b1rms4seq_optFSE_TRvar_FBP(Sample_weight,vector_dTE(jj), TR, ...
                        vector_ETL(ii), vector_flipAngle(gg), sliceThickn, nslices, res,'Fals');
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
    
    % ... 3.4 - Save data ...
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
end


if testSAR == 'loadD'
    
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
        cd(file_path)
    end
end   
toc

%Models_accepted: [TE(ms) #ETL angle TimeScan(min) b1+rms(uT) Trec(ms)]
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
                                            FA_exc_dic, B1, T1maxKnee,  vector_T2(gg),...
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
% %                 fprintf(['      uCRLB ',num2str(uCRLB(ii,jj,gg)),'\n'])
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
SNR_val   = vector_SNR(1);     % Possible SNR values = 1 30 100
T2_val    = vector_T2(2);     % units (ms) - Possible T2 values = 8    16    24    32    40    48
testdTE   = 8;     % units (ms)
testETL   = 10;      % echoes
testFlipA = vector_flipAngle(13);    % units (º) - 115

SNRindx    = find(vector_SNR == SNR_val);
aux_T2indx = find(vector_T2 == T2_val);
T2indx     = (size(Models_accepted,2)+aux_T2indx);
% T2thrs     = 10*T2_val;

% ... 5.1 - Get  Results ...
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
intervPlot_vardT2 = [0 max(max(log(Results_vardT2(:,T2indx:end))))]; % for caxis

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
    % write in excell the parameters    
    B = {'Test','log(T2 max variance)','T2 (ms)','dTE (ms)','ETL','SNR','flipAngle (º)','Trec (s)','TR_acq (s)'};
    sheet = 1+gg;
    xlRange_title = 'A1';
    xlswrite(filename_excell,B,sheet,xlRange_title)
    
    fprintf(['-------------------\n\n'])
    
    [maxValue_vardT2(gg), indxResult_vardT2] = (max( (Results_vardT2(:,T2indx-aux_T2indx+gg)) ) ) ;
    fprintf(['   T2 max variance ' ,num2str(log(maxValue_vardT2(gg))),' for: ', ...
        'T2 = ', num2str(vector_T2(gg)),...        
        ' ms, dTE = ', num2str(Results_vardT2(indxResult_vardT2,1)),...
        ' ms, ETL = ', num2str(Results_vardT2(indxResult_vardT2,2)),...
        ' SNR = ',num2str(SNR_val),...
        ', flipAngle = ', num2str(Results_vardT2(indxResult_vardT2,3)),'º',...
        ', Trec = ', num2str(Results_vardT2(indxResult_vardT2,7)*1e-3),'s',...
        ', TR_acq = ',num2str(Results_vardT2(indxResult_vardT2,8)),'s\n\n'])
    
    results.dTE_vector(gg)       =  Results_vardT2(indxResult_vardT2,1);
    results.ETL_vector(gg)       =  Results_vardT2(indxResult_vardT2,2);
    results.SNR_val_vector(gg)   =  SNR_val;
    results.flipAngle_vector(gg) =  Results_vardT2(indxResult_vardT2,3);
    results.TR_acq_vector(gg)    =  Results_vardT2(indxResult_vardT2,7);
    results.TR_acq_vector(gg)    =  Results_vardT2(indxResult_vardT2,8);
    
    fprintf(['-------------------\n'])

    % write in excell the parameters
    C = {testFSEopt,log(maxValue_vardT2(gg)),vector_T2(gg),Results_vardT2(indxResult_vardT2,1), ...
            Results_vardT2(indxResult_vardT2,2),SNR_val,...
            Results_vardT2(indxResult_vardT2,3),...
            Results_vardT2(indxResult_vardT2,8)};
    xlRange_T2 = ['A',num2str(testFSEopt)];
    xlswrite(filename_excell,C,sheet,xlRange_T2)    
end


%% 7 - Figures (ISMRM23 - Plots)
figure()

for aa= 1:size(vector_T2,2)
    clear aux_results_dTE aux_results_ETL aux_results_flipA a b c uni_dTE_x uni_dTE_y ...
        uni_ETL_x uni_ETL_y uni_flipA_x uni_flipA_y results_dTE results_ETL results_flipA ...
        Results

    % ... 5.1 - Parameters ...
    SNR_val   = results.SNR_val_vector(aa);      % Possible SNR values = 1 30 100
    T2_val    = vector_T2(aa);                   % units (ms) - Possible T2 values = 8    16    24    32    40    48
    testdTE   = results.dTE_vector(aa);          % units (ms)
    testETL   = results.ETL_vector(aa);          % echoes
    testFlipA = results.flipAngle_vector(aa);    % units (º) - 115

    SNRindx    = find(vector_SNR == SNR_val);
    aux_T2indx = find(vector_T2 == T2_val);
    T2indx     = (size(Models_accepted,2)+aux_T2indx);
    % T2thrs     = 10*T2_val;

    % ... 5.1 - Get  Results ...
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
    % plot for fixed dTE
    subplot(3,1,aa)
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
    hold on

end
intervPlotISMRM23 = [0 max(max(log(Results_vardT2(:,:))))]; % for caxis
caxis(intervPlotISMRM23)
