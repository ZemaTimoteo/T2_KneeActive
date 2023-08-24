% Get derivatives

step    = 1e-4;


%% 1 - Derivative from myEPG.m

% --- 1.1 S(M0,T2+step) ---
S                  = All_dict_pp;
[S_T2_plus_step,~] = auxOptimiz_dict_pars_generator_ASM(T1_dic,T2_dic+step,B1_dic,...
                                     dTE,ETL,phase_refoc,phase_exc,FA_exc_dic, FA_refoc_dic,[],methodDic); % 10min

% auxVariab_S_T2 = exp( - Trec/T1); % control for Trecovery & TR variation
% aux2var_S_T2   = (1 - auxVariab_S_T2) / ( 1 - auxVariab_S_T2 * S_T2_plus_step(ETL) );
% S_T2_plus_step = aux2var_S_T2.* S_T2_plus_step;

% --- 1.2  Num derivatives (forward differentiation) ---
dS_dT2_myEPG = (abs(S_T2_plus_step) - abs(S))./(step);

% --- 1.3 Figures ---
% % figure, subplot(121), plot(abs(S),'*--'), hold on, plot(abs(S_T2_plus_step),'ro--')
% % legend('Original Signal','T2 + step')
% % subplot(122), plot(abs(dS_dT2_myEPG),'g+--'), hold on
% % legend('Derivative')


%% 2 - Derivative from Shahin's code

% Get EPG
testS       = testF;
[dF_dT2]    = auxOptimiz_EPG_forward1(alpha0,phase0,'ESP',dTE,'T1',T1_dic,'T2',T2_dic+step); % set target, first simulate for ideal alpha
testdS_dTe2 = dF_dT2(1,3:end);

% --- 1.2  Num derivatives (forward differentiation) ---
dS_dT2_Shahan = (abs(testdS_dTe2) - abs(testS))./(step);