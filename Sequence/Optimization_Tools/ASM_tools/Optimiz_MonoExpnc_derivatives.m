%% FSE optimization with CRLB analysis fitting mono-exponential curves
% With mono-exponential approach - numerical derivatives and analitical derivatives
% TTFernandes - Sept2021

%% 0 - Set matlab paths and toolboxes
file_path = 'D:\Tiago\Trabalho\2021_2025_PhD\Projects\qMRI_Joint\Code\matlabCode\qMRI_tools\Sequence_Optimization';
addpath(genpath(file_path))

%% 1 - Tests
plotTest    = 'True';
saveResults = 'True';

testDict    = 'True';    % Test Dictionary - 'True' OR Load Dictionary - 'False'
testMC      = 'True';    % Test Monte-Carlo for sanity check of CRLB model
T2test      = 'indepnt'; % 'M0depnd' - It is dependent from M0 OR 'indepnt' it's independent from M0
GaussFit    = 'Fals';


%% 3 - Get T2 Signal - Mono-exponencial 

% ... 3.1 - Parameters for generating signals with slr_profile.py and my_epg.py ...
step    = 1e-4;
S0      = 1;   % equal to M0
TE      = [dTE:dTE:dTE*ETL];
% SNR = 30;

% ... 3.2 - Signal S - Mono-Exponencial ...
S     = S0 * exp(-TE/T2);
for ii=1:length(SNR)
    sigma(ii) = (S0/SNR(ii));  % Noise variance for CRLB - Uma diferença
end

% ... 4.3 Analitical derivatives ...
dS_M0_analit = exp(-TE/T2);
dS_T2_analit = S0*((TE)/(T2.^2)).*exp(-TE/T2);  % derivitavite in T2

% % if plotTest == 'True'
% %     figure, plot(abs(S),'*--'), hold on, plot(abs(dS_M0_analit),'ro--'), hold on, plot(abs(dS_T2_analit),'g+--')
% %     hold on
% %     title('Analitical derivatives')
% %     legend('Original Signal','M0 derivative','T2 derivative')
% % end

% ... 4.4 Numerical derivatives ...
S_M0_plus_step = (S0+step) * exp(-TE/T2);     % --- S(M0+step,T2) ---
S_T2_plus_step = S0 * exp(-TE/(T2+T2*step));     % --- S(M0,T2+step) ---

dS_M0_numer      = (abs(S_M0_plus_step) - abs(S))./(step);     % Num derivatives M0 (forward differentiation) ---
dS_T2_numer = (abs(S_T2_plus_step) - abs(S))./(step);     % Num derivatives T2 (forward differentiation) ---


%% 5 - CRLB
% ... 5.1 Jacobian Matrix ...
dM_ana      = [dS_M0_analit(:) dS_T2_analit(:)];     
dM_num      = [dS_M0_numer(:)  dS_T2_numer(:)]; 

% ... 5.2 CRLB ...
tic
for ii=1:length(SNR)    
    FIM_ana       = dM_ana.'   *  (   (1/sigma(ii)^2) * eye( size(dM_ana,1) )  )   *  dM_ana;
    CRLB_ana{ii}  = diag(  sqrt( inv(FIM_ana)) );
    uCRLB_ana(ii) = abs(sqrt(CRLB_ana{ii}(2)))/T2;    % in Zhang et al. 2013 in Proc. EMBS
    
    FIM_num       = dM_num.'   *  (   (1/sigma(ii)^2) * eye( size(dM_num,1) )  )   *  dM_num;
    CRLB_num{ii}  = diag(  sqrt( inv(FIM_num) )  );
    uCRLB_num(ii) = abs(CRLB_num{ii}(2))/ T2;    % in Zhang et al. 2013 in Proc. EMBS
end
toc
% % fprintf(['      uCRLB_ana = ',num2str(uCRLB_ana(1)*sqrt(ETL*dTE)),...
% %     ' / T2 = ',num2str(T2),' / ETL = ',num2str(ETL),' / dTE = ',num2str(dTE),'\n\n'])
% % fprintf(['      uCRLB_num = ',num2str(uCRLB_num(1)*sqrt(ETL*dTE)),...
% %     ' / T2 = ',num2str(T2),' / ETL = ',num2str(ETL),' / dTE = ',num2str(dTE),'\n\n'])

fprintf(['      uCRLB_ana = ',num2str(min(uCRLB_ana)),...
    ' / T2 = ',num2str(T2),' / ETL = ',num2str(ETL),' / dTE = ',num2str(dTE),'\n\n'])
fprintf(['      uCRLB_num = ',num2str(min(uCRLB_num)),...
    ' / T2 = ',num2str(T2),' / ETL = ',num2str(ETL),' / dTE = ',num2str(dTE),'\n\n'])
