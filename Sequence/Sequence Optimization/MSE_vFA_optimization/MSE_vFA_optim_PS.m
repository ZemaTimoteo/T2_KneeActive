%% MSE optimization parameters for variable Flip Angle (vFA) 
% Optimized with Pattern Search (PS) for CRLB with constraint of SNR, SAR and Time FSE for specific T2 values
% by; TTFernandes, Aug - 2024
% Pulseq version: pulseq-1.4.0
% @IST
% Needs Genetic algorithm toolbox from Matlab

clear all
clc
close all

%% Settings
% ... Set Paths ...

% Active path
filePath                 = matlab.desktop.editor.getActiveFilename;
justPath                 = find(filePath == 'M');
DIR_sequenceOptimization = filePath([1:justPath(end)-1]); clear filePath
idx_tb_Path              = find(DIR_sequenceOptimization == 'M');
toolbox_Path             = [DIR_sequenceOptimization([1:idx_tb_Path(end)-1]) 'Toolboxes']; 
cd(DIR_sequenceOptimization)

DIR_data                 = [DIR_sequenceOptimization 'Data'];
DIR_resultsOptimization  = [DIR_data filesep 'PS_sequenceOptimization'];

addpath(genpath(DIR_sequenceOptimization))
addpath(genpath(toolbox_Path))

% installRFtools_mex

%% 0 - Tests 
% ... 0.1 - Select to choose w/ or w/out Constraints ...
testSequence  = 'True';     % Test Optimization OR Read Data (already runned): 'True' or 'Load'

testContrains = 'True';     % Test with Constraints: 'True' or 'Fals'
testIni       = 'True';     % Check initial Conditions: 'True' or 'Fals'
testX0        = 'fixFA';    % 'varFA' OR 'fixFA' OR 'preFA' OR 'conFA'
testDict      = 'JUSTdict'; % SLR profile | For EPG & CRLB - 'SLR_Prof' OR 'JUSTdict'

% ... 0.2 - Set T2 to test ...
test    = 70;
T2_test = [45];  % Short=8(ms), Knee=45(ms), Long=200(ms), GM=83(ms) @3T

%% 1 - Start Code
% ... 1.1 - initialization ...
% TEtest     = 10;    % TE test in (ms)
TEtest     = 10;    % TE test in (ms)
FAini      = 90;     % (º)
FAend      = 60;     % (º)
FixFA      = 160;

% FixFA      = 158;
% vector_ETL = [15:18];    % Set ETL vector - vFA T=45
% vector_ETL = [6:19];    % Set ETL vector
% % vector_ETL = 10;    % Set ETL vector
% % vector_ETL = [6 14];    % Set ETL vector
% vector_ETL = [6:25];    % Set ETL vector
% vector_ETL = [18:22];    % Set ETL vector
% vector_ETL = [6:22];    % Set ETL vector
vector_ETL = [10:25];    % Set ETL vector

% ... 1.2 - Set parameters ...
params               = [];
params.Rayleigh      = 'True';                 % Test with noise 'True' or 'Fals'
params.ParallelImg   = 'LORAK';                % for GRAPPA 'GRAPP' | for LORAKS 'LORAK'
params.RFqual        = 'True';                 % Quality of RF (different slice Thickness in excitation and refocusing - Kumar et al. 2016 JMRI)
params.methodDic     = testDict;               % For EPG & CRLB - 'SLR_Prof' OR 'JUSTdict' OR '90DirSLR' OR '90SLRDir'

params.dir           = DIR_sequenceOptimization;
params.dir_rf        = [DIR_data filesep 'rf_pulses'];
params.dir_data      = DIR_data;

params.T2            = T2_test;                  % T2 (ms) [8 45 200 83]
params               = CF_vFA_PS_parameters(params); % Set parameters for Cost Function (CF)


% ... 1.3 - Cost function ...
fun     = @CF_PS_CRLB_epg_optFSE; % function


%% 2 - Constrains
% SNR > epsolon             - epsolon a pre-defined value (NAO considerar para já)
% SAR < beta                - beta is a pre-defined value
% TE  > 8ms
% TR  > TE * Nx * Nslices
% FA  > 0
% ETL - Another non linear equality constraint could be by using sine or cosine function, e.g sin(k*pi) = 0. This constraint is satisfied only if k is an integer number.

% ... 2.1 - linear constraints
A = [];
b = [];
Aeq = [];
beq = [];

% ... 2.2 - nonlinear constraints
if testContrains == 'True'
    nonlcon = @NonLinearConstrains_PS;
else
    c            = [];
    ceq          = [];
end

%% 3 - Initial value cost function
if testIni == 'True'    
    % 3.1 - Initial Guess ...
    ETLInitest  = vector_ETL(1);    
    if testX0 == 'varFA'
        FA_aux = FAini:(FAend-FAini)/ETLInitest:FAend - (FAend-FAini)/ETLInitest;
    elseif testX0 == 'fixFA'
        FA_aux = ones(1,ETLInitest)*FixFA;
    elseif testX0 == 'preFA'         % results from Oct 24
        FA_aux     = [128.2763      85.79569      81.27881      84.44883      95.77001      113.2271      138.3764      156.9631      161.9245      158.9595       145.513      133.7262      129.1014      134.0628      149.5056      150.5692       148.947];
        TEtest     = 6.5;
        ETLInitest = size(FA_aux,2);
    elseif testX0 == 'conFA'         % results from optimized with constFA
        TEtest     = 7;    
        ETLInitest = 14;
        FA_aux     = ones(1,ETLInitest)*168;
    end
    FA_aux = FA_aux*pi/180; % (rad)    
    x0_aux = [TEtest FA_aux];
    
    % 3.2 -  Initial Value
    Ini_CRLB_aux   = -1 * fun(x0_aux,ETLInitest,params);    % negative because we want to maximize - 'adapted' because is multiplied by negative
    [caux,ceqaux]  = nonlcon(x0_aux,ETLInitest,params);

    % 3.2.2 - Initial for constantFA MSE Siemens 
%     Ini_CRLB_aux   = -1 * fun([13.8 pi*ones(1,5)],5,params);    % negative because we want to maximize - 'adapted' because is multiplied by negative


    % 3.3 - Values
    Ini_AUC          = Constraint_test_AUC(x0_aux,ETLInitest,params);
    [~,Ini_TimeScan] = Constraint_test_Time(x0_aux,ETLInitest,params);
    Ini_SAR          = Constraint_test_SAR(x0_aux,ETLInitest,params); 
    
    % 3.4 - Display
    disp(['Optimized TE               : ' num2str(x0_aux(1)) ])
    disp(['FA vector                  : ' num2str(x0_aux(2:end)*180/pi)])
    disp(['Initial AUC                : ' num2str(Ini_AUC) ])
    disp(['Initial Ini_TimeScan       : ' num2str(Ini_TimeScan) ' (min)'])
    disp(['Initial SAR                : ' num2str(Ini_SAR)])
    disp(['Objective value is         : ' num2str(Ini_CRLB_aux)])
    disp(['Constrains value | diffSAR : ' num2str(caux(1)) ])
    disp(['             ... | diffTscan (s) : ' num2str(caux(2)) ])
    disp(['             ... | diffAUC : ' num2str(caux(3))])

    if testX0 == 'preFA'
        return
    elseif testX0 == 'conFA'
        return
    end
end


%% 4 - Test or Load Sequence Optimization with PS

if testSequence == 'True'
    
    % itterate over different ETL   
    tic
    disp(['Start Pattern Search loop (on T2var|CRLB) over vector ETL'])
    
    parfor kk = 1:size(vector_ETL,2)
        tic
        disp(['ittETL' num2str(vector_ETL(kk))])
        
        
        % 4.1 - Inputs
        
        % ... 4.1.1 - Paramets ...
        ETL = vector_ETL(kk);
        if testX0 == 'varFA'            
            FA_test = FAini:(FAend-FAini)/ETL:FAend - (FAend-FAini)/ETL;
        elseif testX0 == 'fixFA'
            FA_test = [FixFA*ones(1,vector_ETL(kk))];
        end
        FA_test = FA_test*pi/180; % (rad)
        
        % ... 4.1.2 - variable bounds
        lb = [params.constr.betaMin, params.constr.alphaMin*ones(1,ETL)];    % lower bounds - [TE=8(ms),  FA=0º(rad)]
        ub = [params.constr.betaMax, params.constr.alphaMax*ones(1,ETL)];    % upper bounds - [TE=30(ms), FA=180º(rad)]
        
        % 4.2 - Initial Guess ...
        x0{kk}    = [TEtest,FA_test]; % initial guess for: Beta, Alpha - [TE (ms), FA(rad)]        
        if testIni == 'True'    
            Ini_CRLB(kk) = fun(x0{kk},ETL,params); % negative because we want to maximize - 'adapted' because is multiplied by negative
        end
        
        % 4.3 - Options for fmincon
        % ... 4.3.1 - Parameters ...
        rng default % For reproducibility
        
        % ... 4.3.2 - options for the Augmented Lagrangian Genetic Algorithm (ALGA) ...
        options = optimoptions('patternsearch','Display','iter','PlotFcn',@psplotbestf);

            % Mesh
%         options = optimoptions(options,'InitialMeshSize',1/2);  % Initial Mesh Size
        options = optimoptions(options,'InitialMeshSize',1/2);  % Initial Mesh Size
        options = optimoptions(options,'ScaleMesh',false );
        options = optimoptions(options,'AccelerateMesh',true);
%         options.MeshTolerance = 1e-7;       
        options.MeshTolerance = 1e-7;       
        
        % 4.4 - Call Pattern Search
        
        if testContrains == 'True' % With Constraints
            [vector_results{kk},fval(kk),exitflag(kk),outputStruct{kk}] = ...
                patternsearch(@(x) fun(x,ETL,params),x0{kk},A,b,Aeq,beq,lb,ub,@(x) nonlcon(x,ETL,params),options);
        else % Without constraints
            [vector_results{kk},fval(kk),exitflag(kk),outputStruct{kk}] = ...
                patternsearch(@(x) fun(x,ETL,params),x0{kk},A,b,Aeq,beq,lb,ub,[],options); % optimize with Pattern Search
        end
        
        time(kk) = toc;
                
        % 4.5 - Get AUC & TR values
        AUC(kk)             = Constraint_test_AUC(vector_results{kk},ETL,params);
        [TR(kk),Tscan(kk)]  = Constraint_test_Time(vector_results{kk},ETL,params);
        SAR(kk)             = Constraint_test_SAR(vector_results{kk},ETL,params);
        
        % 4.6 - Displays & Print Solution
        % ... 4.6.1 - Display ...
        disp(['................'])
        disp(['itt' num2str(kk)])
        disp(['Time: ' num2str(time(kk)) ' (s)'])
        disp(['................'])        
        disp(['Objective value (T2var|CRLB) is: ' num2str(-fval(kk))])
        
        % ... 4.6.2 - Print solution ...
        disp(['................'])
        disp('Solution')
        disp(['T2  = ' num2str(params.T2),   ' (ms)'])
        disp(['TE  = ' num2str(vector_results{kk}(1)), ' (ms)'])
        disp(['ETL = ' num2str(ETL)])
        disp(['FA  = ' num2str(vector_results{kk}(2:end)*180/pi), ' (º)'])
        disp([''])        
        disp(['...Constraints...'])
        disp(['AUC     = ' num2str(AUC(kk)),   ' (a.u.) | minVal = ' , num2str(params.constr.minAUC) ])
        disp(['B1+rms  = ' num2str(SAR(kk)),   ' (uT)   | maxVal = ' , num2str(params.constr.maxB1_rms) ])                
        disp(['TR      = ' num2str(TR(kk)),    ' (s)'  ])                
        disp(['Tscan   = ' num2str(Tscan(kk)), ' (min)  | maxVal = '   , num2str(params.constr.T_acq)])                
        disp(['................'])
        
    end
    
    
    %  4.7 Save Results
    fval                   = - fval;          % CRLB - make it the maximum | multiply by (-) to get maximum
    results.vector_results = vector_results;
    results.fval           = fval;            
    results.exitflag       = exitflag;
    results.outputStruct   = outputStruct;
    results.time           = time;
    results.x0             = x0;
    results.params         = params;
    results.AUC            = AUC;
    results.vector_ETL     = vector_ETL;
    results.TR             = TR;
    results.SAR            = SAR;
    results.Tscan          = Tscan;
    
    cd(DIR_resultsOptimization)
    save(['PS_optimizationSequence_CRLB_test'  num2str(test) '_T2_' num2str(params.T2) '.mat'],'results')
    
    
elseif testSequence == 'Load'     %  4.8 Load Results

    cd(DIR_resultsOptimization)    
    load(['PS_optimizationSequence_CRLB_test'  num2str(test) '_T2_' num2str(params.T2) '.mat'])
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
    
end

beep


% test for t2=45ms
if T2_test==45 && test == 9
    idx = 1:9;
    fval           = fval(idx);
    time           = time(idx);
    x0             = x0(idx);
    vector_results = vector_results(idx);
    AUC            = AUC(idx);
    vector_ETL     = vector_ETL(idx);
    TR             = TR(idx);
    SAR            = SAR(idx);
    Tscan          = Tscan(idx);
end
%% 5 - Find best Results

% 5.1 - Organize Results
optimizedFA = NaN(size(vector_ETL,2),(vector_ETL(end)));
for kk = 1:size(vector_ETL,2)   
    optimizedTE(kk)                  = real(vector_results{kk}(1));
    optimizedFA(kk,1:vector_ETL(kk)) = real(vector_results{kk}(2:end)*180/pi);
    
    % %     fprintf([' ----------------Best Solution per K-------------\n']);
    % %     disp(['Value for xTest of:         - TE  = ', num2str(optimizedTE(kk)), ' (ms)']);
    % %     disp(['                            - ETL = ', num2str(vector_ETL(kk))]);
    % %     disp(['                            - FA  = ', num2str(optimizedFA(kk,1:vector_ETL(kk))), ' (º)']);
    % %     disp(['Cost Function (vardT2) value is   : ', num2str(fval(kk))]);
    % %     fprintf([' --------------------------------------------\n']);
end

% 5.2 - Find best Solution
[~,idx] = max(fval);    
fprintf([' ----------------Best Solution----------------\n']);
disp(['Value for xTest of:         - TE    = ', num2str(optimizedTE(idx)), ' (ms)']);
disp(['                            - ETL   = ', num2str(vector_ETL(idx))]);
disp(['                            - FA    = ', num2str(optimizedFA(idx,1:vector_ETL(idx))), ' (º)']);
disp(['                            - TR    = ', num2str(TR(idx)), ' (s)']);
disp(['                            - Tscan = ', num2str(Tscan(idx)), ' (min)']);
disp(['Cost Function (T2var|CRLB) is       : ', num2str(fval(idx))]);
disp(['AUC value                           : ', num2str(AUC(idx))]);
disp(['SAR                                 : ', num2str(SAR(idx))]);
fprintf([' --------------------------------------------\n']);

beep

%% 6 - Generate EPGs
figure()
for ii=1:size(vector_ETL,2)
    ETL            = vector_ETL(ii);
    TE             = vector_results{ii}(1);
    timeVector{ii} = [TE:TE:TE*ETL];    
    data           = CF_epg_optFSE_getCRLB(vector_results{ii},ETL,params);
    signal{ii}     = data.signal;
    
    if ii == idx % bigger line for best fit
        plot(timeVector{ii} ,signal{ii},'LineWidth',5)
    else
        plot(timeVector{ii} ,signal{ii})
    end
    hold on
end
xlabel('Time (ms)')
xline(T2_test*1.6)
title(['EPG of each ETL studied T2=' num2str(params.T2)])
legend([num2str(ETL)])
% % legend(['ETL15'],['ETL16'],['ETL17'],['ETL18'], ...
% %             ['Measure of AUC for 150% of T2=45 (ms)'])
legend(['ETL6'],['ETL7'],['ETL8'],['ETL9'],['ETL10'],['ETL11'],['ETL12'],['ETL13'],['ETL14'],['ETL15'],['ETL16'],['ETL17'],['ETL18'],['ETL19'], ...
            ['Measure of AUC for 150% of T2=45 (ms)'])
% % legend(['ETL8'],['ETL9'],['ETL10'],['ETL11'],['ETL12'],['ETL13'],['ETL14'],['ETL15'],['ETL16'],['ETL17'],['ETL18'],['ETL19'],['ETL20'],['ETL21'],['ETL22'], ...
% %             ['Measure of AUC for 150% of T2=45 (ms)'])
% % legend(['3'],['4'],['5'],['6'],['7'],['8'],['9'],['10'],['11'],['12'],['13'],['14'],['15'],['16'],['17'],['18'],['19'],['20'],['21'],['22'],[ ...
% %             '23'],['24'],['25'],['26'],['27'],['28'],['29'],['30'])


% % T2_83_vFA = signal{idx};
% % T2_45_vFA = signal{idx};
% % T2_8_vFA = signal{idx};

%% 7 - Displays & Print Solution

fprintf(['\n \n 7 - Displays & Print Solution for T2= ' num2str(T2_test) 'ms  \n \n'])

% options

% ... 7.1 - Initial Display ...
% % if testIni == 'True'
% %     disp([' ---------Optimz. w/ Pattesearch ----- '])
% %     disp(['ETL: ' num2str(vector_ETL(idx))])
% %     disp(['Initial - TE: ' num2str(x0{idx}) ' (ms); & FA vector: ' num2str(x0{idx}(2:end)*180/pi) ' (º)'])
% %     disp(['Objective value (T2var|CRLB) is: ' num2str(Ini_AUC(kk))])
% % end
disp([' ------------------------------------- '])    

% ... 7.2 - Optimized Display ...
disp(['ETL: ' num2str(vector_ETL(idx))])
disp(['Optimized - TE: ' num2str(optimizedTE(idx)) '; & FA vector: ' num2str(optimizedFA(idx,1:vector_ETL(idx)))])
disp(['Objective value (T2var|CRLB) is: ' num2str(fval(idx))])
disp(['AUC (T2var|CRLB) is: ' num2str(AUC(idx))])
disp(['SAR_B1rms+ (uT) is: ' num2str(SAR(idx))])
disp(['Tscan (min) is: ' num2str(Tscan(idx))])
disp(['TR (s) is: ' num2str(TR(idx))])
disp(['Running Time: ' num2str(time(idx))])

disp([' ------------------------------------- '])

beep
%% 8 - Constraints

% .. 8.1 - Initial Values ...
fprintf('\n \n 8 - Get Constraints \n \n')
disp(['           CRLB val | value' ...
    ' | Constraints (Max Values) | SAR: ' num2str(params.constr.maxB1_rms) ...
                          '(uT) | TimeScan: ' num2str(params.constr.T_acq) ...
                         '(min) | AUC: ' num2str(params.constr.minAUC)])

% .. 8.2 - Initial Differences ...                         
if testIni == 'True'
    [Ini_c,~] = nonlcon(x0{idx},vector_ETL(idx),params);                 % Values of contrains for x0
    disp(['Ini        ' ...
           'CRLB val | ' num2str(params.constr.minVarT2)...
            '| Constraints | Diff SAR: ' num2str(Ini_c(1)) ...
                        '(uT) | Diff TimeScan: ' num2str(Ini_c(2)) ...
                        '(s)  | Diff AUC:' num2str(Ini_c(3)) ...
                            ' | (All need to be -)'])
end

% .. 7.3 - Optimized Differences ...                         
[End_c,~] = nonlcon(vector_results{idx},vector_ETL(idx),params);      % NEED to BE NEGATIVE
disp(['Optimiz    ' ...
    'CRLB val | ' num2str(fval(idx))...
    '| Constraints | Diff SAR: ' num2str(End_c(1)) ...
                    '(uT) | Diff TimeScan: ' num2str(End_c(2)) ...
                    '(s)  | Diff AUC:' num2str(End_c(3)) ...
                        ' | (All need to be -)'])
                 
beep                    
%% 9 - Final Figures

fprintf('\n \n 9 - Get Figures  \n \n')
 
% 9.1 - TE Figure
figure()
subplot(1,3,1)
p = plot(optimizedTE);
p.LineWidth = 2;
hold on
title(['Optimized TE - T2=' num2str(params.T2) ' ms'])
xlabel('# ETLs','fontsize',15)
ylabel('TE (ms)','fontsize',15)
set(gca,'FontSize',10)
xticks(1:1:size(vector_ETL,2)*2)
xticklabels(vector_ETL)

% 9.2 - FA Figure
plotOptimFA = zeros(size(optimizedFA,1)+1,size(optimizedFA,2));
plotOptimFA(1:size(optimizedFA,1),:) = optimizedFA;
subplot(1,3,2)
pcolor(plotOptimFA)
hold on
title('Optimized FA (º)')
xlabel('# FA','fontsize',15)
ylabel('ETL','fontsize',15)
set(gca,'FontSize',10)
xticks(1:1:size(vector_ETL,2)*4)
yticks(1:1:size(vector_ETL,2)*2+2)
xticklabels(1:max(vector_ETL))
yticklabels(vector_ETL)
caxis([min(min(optimizedFA))-5 max(max(optimizedFA))+5])
colorbar()
colormap(hot)

% 9.3 - Cost Function Figure
subplot(1,3,3)
p = plot(fval,'r');
p.LineWidth = 2;
hold on
title(['Cost Function (T2var|CRLB) | T2= ' num2str(params.T2) ' ms'])
xlabel('# ETLs','fontsize',15)
set(gca,'FontSize',10)
xticks(1:1:size(vector_ETL,2)*2)
xticklabels(vector_ETL)
ylabel('dE_dT2','fontsize',15)

%% Plot For ISMRM_2024
% 9.1 - TE Figure
figure()
subplot(1,2,1)

yyaxis left
p = plot(vector_ETL,optimizedTE);
p.LineWidth = 2;
hold on
ylabel('TE (ms)','fontsize',25)
set(gca,'FontSize',20)

yyaxis right
p = plot(vector_ETL,fval,'r');
p.LineWidth = 2;
hold on
xlabel('# ETLs','fontsize',25)
set(gca,'FontSize',20)
ylabel('T2var | CRLB','fontsize',25)

title(['vFA | Optimized TE and Cost Function for T2= ' num2str(params.T2) ' ms'])
legend(['TE (ms)'],['Cost Function (T2var | CRLB)'])

% 9.2 - FA Figure
plotOptimFA = zeros(size(optimizedFA,1)+1,size(optimizedFA,2));
plotOptimFA(1:size(optimizedFA,1),:) = optimizedFA;
subplot(1,2,2)
pcolor(plotOptimFA)
hold on
title(['vFA | Optimized FA (º)' ' for T2= ' num2str(params.T2) ' ms'])
xlabel('# FA','fontsize',25)
ylabel('ETL','fontsize',25)
set(gca,'FontSize',20)
xticks(1:1:size(vector_ETL,2)*4)
yticks(1:1:size(vector_ETL,2)*2+2)
xticklabels(1:max(vector_ETL))
yticklabels(vector_ETL)
caxis([min(min(optimizedFA))-5 max(max(optimizedFA))+5])
colorbar()


