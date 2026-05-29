%% Test CRLB Values for PS
% by; TTFernandes, Oct - 2023
% @IST

clear all 
clc

%% 0 - Settings
% ... 0.1 - Set Paths ...
filePath                 = matlab.desktop.editor.getActiveFilename;
justPath                 = find(filePath == 't');
DIR_sequenceOptimization = filePath([1:justPath(end-1)-1]); clear filePath
idx_tb_Path              = find(DIR_sequenceOptimization == 'M');
toolbox_Path             = [DIR_sequenceOptimization([1:idx_tb_Path(end)-1]) 'Toolboxes']; 
cd(DIR_sequenceOptimization)

addpath(genpath(DIR_sequenceOptimization))
addpath(genpath(toolbox_Path))

DIR_resultsOptimization = [DIR_sequenceOptimization 'Data' filesep 'PS_sequenceOptimization'];

%% Get parameters
% Test Parameters
test                 = 31;
T2_test              = 45;

params               = [];
params.Rayleigh      = 'True';                 % Test with noise 'True' or 'Fals'
params.ParallelImg   = 'GRAPP';                % for GRAPPA 'GRAPP' | for LORAKS 'LORAK'
params.RFqual        = 'True';                 % Quality of RF (different slice Thickness in excitation and refocusing - Kumar et al. 2016 JMRI)
testDict             = 'SLR_Prof';             % SLR profile | For EPG & CRLB - 'SLR_Prof' OR 'JUSTdict'
params.methodDic     = testDict;               % For EPG & CRLB - 'SLR_Prof' Test=(test31 T2=45 & test33 - T2=8) OR 'JUSTdict' (test30 - T2=45 & test32 - T2=8)
params.T2            = T2_test;                % T2 (ms) [8 45 200 83]
params.dir_rf        = [DIR_sequenceOptimization 'Data' filesep 'rf_pulses'];
params.dir           = DIR_sequenceOptimization;

% params          = CF_Mycode_parameters_NLO(params); % Set parameters for Cost Function (CF)
params          = CF_PS_parameters(params); % Set parameters for Cost Function (CF)

display(['Dictionary tested with ' testDict])

%% Load - data from test:
loadName = ['PS_optimizationSequence_CRLB_test' num2str(test) '_T2_' num2str(T2_test) '.mat'];
cd(DIR_resultsOptimization)
load(loadName)

vector_results = results.vector_results;
fval           = results.fval;
exitflag       = results.exitflag;
outputStruct   = results.outputStruct;
time           = results.time;
x0             = results.x0;
params         = results.params;
AUC            = results.AUC;
vector_ETL     = results.vector_ETL;
TR             = results.TR;
SAR            = results.SAR;
Tscan          = results.Tscan;

%% Find best Solution
% Organize Results
optimizedFA = NaN(size(vector_ETL,2),(vector_ETL(end)));
for kk = 1:size(vector_ETL,2)
    optimizedTE(kk)                  = real(vector_results{kk}(1));
    optimizedFA(kk,1:vector_ETL(kk)) = real(vector_results{kk}(2:end)*180/pi);
end

% Get Result
[~,idx] = max(fval);
fprintf([' ----------------Best Solution----------------\n']);
disp(['Value for xTest of:         - TE    = ', num2str(optimizedTE(idx)), ' (ms)']);
disp(['                            - ETL   = ', num2str(vector_ETL(idx))]);
disp(['                            - FA    = ', num2str(optimizedFA(idx,1:vector_ETL(idx))), ' (º)']);
disp(['                            - TR    = ', num2str(TR(idx)), ' (s)']);
disp(['                            - Tscan = ', num2str(Tscan(idx)), ' (min)']);
disp(['Cost Function (T2var|CRLB) is       : ', num2str(fval(idx))]);
disp(['AUC value                           : ', num2str(AUC(idx))]);
fprintf([' --------------------------------------------\n']);

beep

%% Test CRLB
for ii=1:size(results.vector_results,2)
    ii
    xtest                   = results.vector_results{1,ii};
    ETLtest(ii)             = size(results.vector_results{1,ii},2)-1;
    CRLBtest(ii)            = - CF_PS_CRLB_epg_optFSE(xtest,ETLtest(ii),params);  % check sensitivity to estimate T2 values    
    [TRtest(ii), Tscan(ii)] = Constraint_test_Time(xtest,ETLtest(ii),params);     % check TR and Time of scan       
    AUCtest(ii)             = Constraint_test_AUC(xtest,ETLtest(ii),params);      % check Area under the curve  
    SARtest(ii)             = Constraint_test_SAR(xtest,ETLtest(ii),params);      % check SAR       
end

%% Best Result
[bestCRLB, idxCRLB] = max(CRLBtest(:));
aux_vector          = results.vector_results{1,idxCRLB};
TEbest              = aux_vector(1);
ETLbest             = ETLtest(idxCRLB);
vFAbest             = aux_vector(2:end)*180/pi;
TRbest              = TRtest(idxCRLB);
Tscanbest           = Tscan(idxCRLB);
AUCbest             = AUCtest(idxCRLB);
SARbest             = SARtest(idxCRLB);

%% Print Best results
fprintf([' ----------------Best Solution----------------\n']);
disp(['Best Parameters for T2=' num2str(T2_test) '(ms) are:   - TE    = ', num2str(TEbest), ' (ms)']);
disp(['                                     - ETL   = ', num2str(ETLbest)]);
disp(['                                     - FA    = ', num2str(vFAbest), ' (º)']);
disp(['                                     - TR    = ', num2str(TRbest), ' (s)']);
disp(['                                     - Tscan = ', num2str(Tscanbest), ' (min)']);
disp(['Cost Function (T2var|CRLB) value is  : ', num2str(bestCRLB)]);
disp(['AUC value                            : ', num2str(AUCbest)]);
disp(['SAR (uT)                             : ', num2str(SARbest)]);
fprintf([' --------------------------------------------\n']);

%% Figuresestamos 
figure()
plot(ETLtest,CRLBtest)
hold on
title(['test ' num2str(test) ' | CRLB evolution'])

%% EPG
figure()
xBest      = results.vector_results{1,idxCRLB};
data       = CF_epg_optFSE_getSignalEPG(xBest,ETLbest,params);
signal = data.signal;

plot(signal,'linewidth',2)
hold on

title(['EPG of each T2 optz. for variable Flip Angles'])
% % legend(['Short | 8    (ms)'],['Knee | 45   (ms)'],['GM | 200 (ms)'],['GM   | 83   (ms)'])
legend(['Knee | 45   (ms)'])
