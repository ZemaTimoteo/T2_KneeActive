%% MSE optimization with CRLB analysis tested with Monte-Carlo Simulations fitting EPG curves
% Dictionaries previously generated
% TTFernandes - Jan2025


% =============================
% Toolboxes:
%      - sar4seq - in: https://github.com/imr-framework/sar4seq/tree/master/sar4seq
% =============================
%% 0 - Set matlab paths and toolboxes
clear, clc
close all


%% 0 - Settings
% ... 0.1 - Set Paths ...
filePath                 = matlab.desktop.editor.getActiveFilename;
justPath                 = find(filePath == '/');
pathSeqOpt               = filePath([1:justPath(end)]);  clear filePath
DIR_sequenceOptimization = [pathSeqOpt 'MSE_vFA_optimization'];
idx_tb_Path              = find(DIR_sequenceOptimization == 'M');
toolbox_Path             = [pathSeqOpt 'Toolboxes']; 
cd(DIR_sequenceOptimization)

addpath(genpath(pathSeqOpt))
addpath(genpath(DIR_sequenceOptimization))
addpath(genpath(toolbox_Path))

%% 1 - Get parameters
% Input Parameters
B1_scanner = 1; % Teste de 11/9

T1_test = 1000;
T2_test = 45;
T2_dic  = [8:80];

params          = [];
params.Rayleigh      = 'True';                 % Test with noise 'True' or 'Fals'
params.ParallelImg   = 'GRAPP';                % for GRAPPA 'GRAPP' | for LORAKS 'LORAK'
params.RFqual        = 'True';                 % Quality of RF (different slice Thickness in excitation and refocusing - Kumar et al. 2016 JMRI)
testDict             = 'SLR_Prof';             % RF profile | 'SLR_Prof' OR 'JUSTdict'
testPlot             = 'Fals';
testGaussFit         = 'Fals';
params.methodDic     = testDict;               % For EPG & CRLB - 'SLR_Prof' OR 'JUSTdict' OR '90DirSLR' OR '90SLRDir'
params.T2            = T2_test;                % T2 (ms) [8 45 200 83]
params.dir_rf        = [DIR_sequenceOptimization filesep 'Data' filesep 'rf_pulses'];
params.dir           = DIR_sequenceOptimization;

% params          = CF_Mycode_parameters_NLO(params); % Set parameters for Cost Function (CF)
params          = CF_PS_parameters(params); % Set parameters for Cost Function (CF)

display(['Dictionary tested with ' testDict])

%% 2 - Test cFA CRLB and AUC - from ISMRM 23 abstract TTFernandes - CRLB minimization for optimal Multi Spin-Echo T2-mapping of the cartilage using dictionary-based estimation
% parameters
% [ Optimz for T2=8ms | T2=45 | T2=83 | Normal MSE params | Study in HSa ]
if testDict == 'JUSTdict'
    T1_ini  = [ones(1,2)*1000   1000    1000  ];
    T2_ini  = [   8      45       45      45  ];
    TE_ini  = [   8.35    6.5     10      13.8];
    ETL_ini = [  12      18        7       5  ];
    FA_ini  = [ 125     123      180     180  ]*B1_scanner;
    TR_ini  = [2420    3673     2700    1910  ];
elseif testDict == 'SLR_Prof'
    T1_ini  = [T1_test     T1_test   T1_test  ];
    T2_ini  = [  45       45     45  ];
    TE_ini  = [   6.5     10     13.8];
    ETL_ini = [  18       11      5  ];
    FA_ini  = [ 124      180    180  ]*B1_scanner;
    TR_ini  = [3161     5000   1910  ];
end

AUCcutTime = min(ETL_ini.*TE_ini);

display(['T2 values used for test (ms)           : ' num2str(T2_test)])
%% 3 - Get initial parameters
CRLBini = zeros(1,size(T1_ini,2));
for jj=1:size(T1_ini,2)
    x_ini_cFA{jj}    = [TE_ini(jj) FA_ini(jj)*ones(1,ETL_ini(jj))*pi/180];      
end

%% 4 - EPG & Dictionary 
% ... 4.1 - Parameters for generating signals with slr_profile.py and my_epg.py ...
T1_dic        = T1_test;                    % ms
B1_dic        = 0.6:0.1:1.4;                % B1 value ( 0.7:0.01:0.75; OR 0.6:0.01:1.4; )
phase_exc     = pi/2;                       % Phase Exciatation - in (rad) - along y
FA_exc_dic    = pi/2;                       % Flip-angle Exciatation - in (rad) - along y
T2_numTest    = size(T2_dic,2);
jj   = 2;                                   % T2=45
Dict = zeros(T2_numTest,ETL_ini(jj));

if testPlot == 'True'
    figure(33)
    xline(AUCcutTime,'LineStyle','--','LineWidth',2,'Color','k')
    hold on
    colorpos = [0 0.4470 0.7410;0.8500 0.3250 0.0980;...
        0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;...
        0.4660 0.6740 0.1880;0.3010 0.7450 0.9330;...
        0.6350 0.0780 0.1840];
end

% 4.2 - Get Dictionary - B1 TODO
for ii = 1:T2_numTest
%     for jj=1:size(T1_ini,2)
    timeVector{jj} = [TE_ini(jj):TE_ini(jj):TE_ini(jj)*ETL_ini(jj)];
    params.T1      = T1_ini(jj);
    params.T2      = T2_dic(ii);
    params.TRini   = TR_ini(jj);

    data        = CF_epg_optFSE_getSignalEPG(x_ini_cFA{jj},ETL_ini(jj),params);
    Dict(ii,:)  = data.signal';

    % ... 4.3 - Figure
    if testPlot == 'True'
        plot(timeVector{jj},signal_CFA{ii,jj}','linewidth',2,'Color',colorpos(jj,:))
        xlabel('Time (ms)')
        hold on
    end
%     end

    fprintf(['Tests for dictionary: ' num2str(ii) ' / ' num2str(T2_numTest) ' \n'])

end

% ... 4.3 - Figure
if testPlot == 'True'
    ylim([0 16])
    title(['EPG of each T2 optz. for constant Flip Angles'])
    % % legend(['Short | 8    (ms)'],['Knee | 45   (ms)'],['GM | 200 (ms)'],['GM   | 83   (ms)'])
    legend(['AUC cut line'],['cFA'],['SE'],['Vender'])
end   

%... 4.4 - Signal
params.T1      = T1_dic;
params.T2      = T2_test;
params.TRini   = TR_ini(1);
data        = CF_epg_optFSE_getSignalEPG(x_ini_cFA{2},ETL_ini(2),params);
Signal      = data.signal;

% ... 4.3 - Figure
if testPlot == 'True'
    figure(34)    
    plot(timeVector{2},Signal','linewidth',2,'Color',colorpos(jj,:))
    xlabel('Time (ms)')
    hold on
end

fprintf('\n\n 4 - Sucessfully finished - Get Dictionary\n\n')


%% 5 - Monte Carlo
% Compare CRLB with standard deviation of Monte Carlo distributions

% ... 5.1 - Parameters for Testing MC...
nreps        = 1e4;               % Number of repetitions
sig          = Signal;            % Signal from EPG
All_dict     = Dict;
All_pars     = T2_dic;
T2_est       = zeros(nreps,1);    % Vectors to store T2 estimates
vector_SNR   = [1 2:2:600];
sigma_vector = zeros(size(vector_SNR));

% ... 5.2 - Normalize dictionary ...
Dict_norm = zeros(size(All_dict,1),size(All_dict,2));
for i=1:size(All_dict,1)
    Dict_norm(i,:) = All_dict(i,:)./norm(All_dict(i,:));
end

% ... 5.3 - Obtain Monte-Carlo (MC) uncertainty ...

for gg=1:size(T2_test,2)
    for ii =1:size(vector_SNR,2) % 67s
        % Parameters
        sigma_vector(ii) = (1/vector_SNR(ii));                   % Noise variance for CRLB        

        % Get CRLB for specific SNR
        params.SNR    = vector_SNR(ii);
        params.sigma3 = sigma_vector(ii);                        
        CRLB_num(ii)  = CF_PS_CRLB_epg_optFSE(x_ini_cFA{2},ETL_ini(2),params);

        % Get MC for the model
% %         fprintf(['SNR test: ', num2str(vector_SNR(ii)),'\n\n'])
        [uMC_fit(ii,gg), uMC_std(ii,gg)] = MC_epg_optMSE(nreps, sig, Dict_norm, All_pars, ...
            sigma_vector(ii), T2_dic, T2_test(gg), CRLB_num(ii), ...
            testGaussFit,testPlot);
        
        fprintf(['MC | T2 values Test: ' num2str(gg) ' / ' num2str(size(T2_test,2)) '| SNR Test: ' num2str(ii) ' / ' num2str(size(vector_SNR,2)) ' \n'])        
    end
end

% ... 5.4 - Figures of MC ...

figure()
plot(vector_SNR,(uCRLB(:,:,3)),'g')
hold on
plot(vector_SNR,(uMC_fit(:,3)),'r')
hold on
plot(vector_SNR,(uMC_std(:,3)),'b')
title('EPG - Uncertainty Comparison')
legend('uncertainty CRLB','uncertainty MC fit','uncertainty MC std')
xlabel('SNR')
ylabel('uncertainty')



fprintf('\n\n 5 - Sucessfully finished - Monte-Carlo estimation\n\n')


