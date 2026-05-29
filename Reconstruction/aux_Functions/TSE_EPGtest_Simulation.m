% Tests T2 MSE for debug of python EPG
% TTFernandes Dec 2023

%% 1 - Set matlab paths and toolboxes
clear all
clc

addpath(genpath('D:\Tiago\Trabalho\2021_2025_PhD\Tutoria\2021-2022_3P_Estagio\Code\matlab\Dictionary_aux_Functions'));
addpath(genpath('D:\Tiago\Trabalho\2021_2025_PhD\Projects\qMRI_Joint\Code\matlabCode\MSE_recon\aux_Functions'));
addpath(genpath('D:\Tiago\Trabalho\2021_2025_PhD\Projects\qMRI_Joint\Code\matlabCode\qMRI_tools\Sequence_Optimization'))

%% 2 - Options
methodDic     = 'DICT_SLR'; % 'DICT_SLR' OR 'JUST_DIC'
dir_Results   = 'D:/Tiago/Trabalho/2021_2025_PhD/Projects/qMRI_Joint/Data/Data_EPG';
dict_dir      = "D:/Tiago/Trabalho/2021_2025_PhD/Projects/qMRI_Joint/Data/Dictionaries";


%% 3 - Get parameters and Signal
% 3.0 - initialization ...
params               = [];
params.dir_rf        = 'D:\Tiago\Trabalho\2021_2025_PhD\Projects\qMRI_Joint\Code\matlabCode\qMRI_tools\Sequence_Optimization\Data\rf_pulses';
params.Rayleigh      = 'True';                   % Test with noise 'True' or 'Fals'
params.ParallelImg   = 'GRAPP';                  % for GRAPPA 'GRAPP' | for LORAKS 'LORAK'
params.methodDic     = 'SLR_Prof';               % For CRLB - 'SLR_Prof' OR 'JUSTdict' OR '90DirSLR' OR '90SLRDir'

dictTest             = 'Fals';
tempMatch            = 'True';                   % 'True' - Get ind_parameters  OR 'Fals' - load ind_parameters

% ... part 3.1 - parameters ....
SNRvalue  = 40;
TE        = 10;                  % Echo Time in (ms)
TR        = 2.500;               % (s)
nslices   = 1;
st        = 3;                   % (mm)
ETL       = 10;                  % Echo Train Length
gamma     = 42.576e6;
Nx        = 256;
Ny        = Nx;
FOV       = 256;
RFexc     = 90;                  % degrees
RFrefoc   = 160;                 % degrees
nc        = 1;                   % Number of Coils for Coil compressing

    % ... 3.1.1 - set Angles ...
FA_exc    = RFexc*pi/180;                          % Flip-angle Exciatation - in (rad) - along y
phase_exc = pi/2;                                  % Phase Exciatation - in (rad) - along y

FA_refoc  = ones(1,ETL) * RFrefoc;                     % (degrees)
phase_refoc = exp(zeros(1,ETL)/(180/pi)*1j);   % Phase Refocosing -  ou zeros(1,ETL);


% ... part 3.2 - Get Sginal ....
T1_signal = 600;  % ms Value from NIST Phantom
T2_signal = 56;  % ms
B1_signal = 1.2;

% ... part 3.3 - get params
params.T2            = 45;                  % T2 (ms) [8 45 200 83]
params               = CF_PS_parameters(params); % Set parameters for Cost Function (CF)
params.T2            = T2_signal;                  % T2 (ms) [8 45 200 83]

[Signal_Evol,pars] = dict_pars_generator_EPGtest_Simulation(T1_signal,T2_signal,B1_signal, ...
                                TE,ETL,phase_refoc,phase_exc,FA_exc, FA_refoc, params, methodDic);

testSignal = Signal_Evol;

% ... part 3.3 - plot ...
figure();
plot(abs(testSignal));

fprintf('\n\n 3 - Sucessfully finished - Define Parametres \n\n')




%% 4 - Build dictionary - bi-exponencial T2 estimation

% ESP - echo Spacing
T1 = T1_signal;                   %ms
T2 = 1:1:100;               %ms
B1 = 0.6:0.01:1.4;

if dictTest == 'True'
    % ESP - Spacing Time between Echoes - in (ms)
    [dict_knee_shortTE,pars] = dict_pars_generator_EPGtest_Simulation(T1,T2,B1,...
                            ESP,nTE,phase_refoc,phase_exc,FA_exc, FA_refoc, params, methodDic); % 10min

    col_T2 = pars(:,1);     % indexes of T2
    col_B1 = pars(:,3);     % indexes of B1
    col_B1 = col_B1*100;    

    % normalize dictionary
    for ii=1:size(dict_knee_shortTE,2)
        Dict_knee_shortTE_norm(:,ii) = dict_knee_shortTE(:,ii)./norm(dict_knee_shortTE(:,ii));
    end

    % Plots
% %         figure()
% %         plot(abs(Dict_knee_shortTE_norm))

    % Save Dictionary
    cd(dict_dir)
    if params.methodDic == 'SLR_Prof'               
        nameMatrix = ['testDICT_SLR_MSE_CRLB_py_TE_', num2str(TE), 'ms_ETL_', num2str(ETL), '_fA_', num2str(RFrefoc)];
    else
        nameMatrix = ['testDICT_dirac_MSE_CRLB_py_TE_', num2str(TE), 'ms_ETL_', num2str(ETL), '_fA_', num2str(RFrefoc)];
    end        
    save([nameMatrix,'.mat'],'Dict_knee_shortTE_norm','pars')

else
    cd(dict_dir)
    if params.methodDic == 'SLR_Prof'               
        nameMatrix = ['testDICT_SLR_MSE_CRLB_py_TE_', num2str(TE), 'ms_ETL_', num2str(ETL), '_fA_', num2str(RFrefoc),'_0'];
    else
        nameMatrix = ['testDICT_dirac_MSE_CRLB_py_TE_', num2str(TE), 'ms_ETL_', num2str(ETL), '_fA_', num2str(RFrefoc),'_0'];
    end
    load(nameMatrix)

    col_B1 = col_B1*100;    % in percentage

end   

fprintf('\n\n 4 - Sucessfully finished - Build dictionary - bi-exponencial T2 estimation\n\n')


%% 5 - T2 Dictionary Matching (selected 5 slices only)
testFORsnr = 1;

if tempMatch== 'True'
    T2_dict = zeros(1, testFORsnr);
    B1_dict = zeros(1, testFORsnr);
    
    for snr=1:testFORsnr
        X = abs(testSignal);
        
        % 5.2.1 - Template match
        ind_param = template_match_EPGtest_Simulation(  dict_Test_norm,  X , col_T2 , col_B1 );
        
        % 5.2.2 - Get T2 and B1 values
        index_Value   = ind_param;
        index_Value   = index_Value.astype(int);
        
        if params.methodDic == 'SLR_Prof'
            T2_dict(0,snr)  = col_T2(index_Value-1)+1;
            B1_dict(0,snr)  = col_B1(index_Value-1)+1;
        else
            T2_dict(0,snr)  = col_T2(index_Value-1)+1;
            B1_dict(0,snr)  = col_B1(index_Value-1)+1;
            
            fprintf('\n\n   -> Itt: ', snr + 1, '/', testFORsnr, 'test')
        end
        
        %# ... 6.3. - Save index ...
        if saveResults
            dir_Results = 'D:/Tiago/Trabalho/2021_2025_PhD/Projects/qMRI_Joint/Data/Data_EPG';
% %             os.chdir(dir_Results)
% %             parametersINDEX    = {'testT2_dict': T2_dict}
% %             parametersINDEX_B1 = {'testB1_dict': B1_dict}

            if SLR_prof  % save results for pre_process
% %                 savemat("testT2_dict_py_SLR.mat", parametersINDEX)
% %                 savemat("testB1_dict_py_SLR.mat", parametersINDEX_B1)
            else
% %                 savemat("testT2_dict_py_Dirac.mat", parametersINDEX)
% %                 savemat("testB1_dict_py_Dirac.mat", parametersINDEX_B1)                    
            end   

        else % load data
            os.chdir(dir_Results)
            if SLR_prof  % save results for pre_process
                Index     = scipy.io.loadmat('testT2_dict_py_SLR');
                Index_B1  = scipy.io.loadmat('testB1_dict_py_SLR');
            else
                Index    = load('testT2_dict_py_Dirac');
                Index_B1 = load('testB1_dict_py_Dirac');
% %                 T2_dict  = Index('T2_dict');
% %                 B1_dict  = Index_B1['B1_dict'];
            end
        end
    end
end


% toc
fprintf('\n\n 5 - Sucessfully finished - Build dictionary - bi-exponencial T2 estimation\n\n')
