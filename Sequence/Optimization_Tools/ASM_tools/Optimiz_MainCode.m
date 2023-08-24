
%% 0 - tests
cd('D:\Tiago\Trabalho\2021_2025_PhD\Projects\qMRI_Joint\Code\matlabCode\qMRI_tools\Optimization_ASM_tests')

methodDic = 'JUST_DIC'; % 'DICT_SLR' or 'JUST_DIC'

%% 1 - Parameters to Optimize Mono-Exponential

% ... 2.1 - Parameters for model optimization ...
T2  = 16;    % T2 values for each we want to optimize FSE sequence

% ... 2.2 - Parameters to model test ...
dTE       = 8;                                  % Echo spacing in ms
ETL       = 4;                                 % Number of Echoes
flipAngle = 180;                                % FLip Angle for Refocusing

% fazer variar o flipangle ao longo da sequencia
SNR    = [15:30:605];                           % variance - % variar SNR de teste
valSNR = SNR(2);
valSigma  = 1/valSNR;                              % value of sigma for CRLB

testPlot = 'Fals';

%% 2 - Parameters for Dictionary
% ... 2.3 - for dictionary ...
T1_dic     = 1000;                  % ms
T2_dic     = T2;        % ms ( 40:1:55; OR 1:1:300; )
B1_dic     = 1;                % B1 value ( 0.7:0.01:0.75; OR 0.6:0.01:1.4; )

phase_exc  = pi/2;                       % Phase Exciatation - in (rad) - along y
FA_exc_dic = pi/2;                       % Flip-angle Exciatation - in (rad) - along y

numflipAngles = flipAngle;                          % Number of Flipe Angles
phase_refoc   = exp(zeros(1,ETL)./(180/pi).*1i);    % phase of refoc pulses (0 = about x axis)
FA_refoc_dic  = ones(1,ETL)*flipAngle;              % Flip-angle Refocosuing - in (degrees) - along y

%% 3 - Test mono_exponential derivatives - numerical & Analytical

% run Mono-Exponential derivatives
Optimiz_MonoExpnc_derivatives

% figures
% % figure, plot(abs(S),'*--'), hold on, plot(abs(dS_M0_analit),'ro--'), hold on, plot(abs(dS_T2_analit),'g+--')
% % title('Numerical derivatives')
% % legend('Original Signal','M0 derivative','T2 derivative')
% % 
    % Plots of derivatives
% % figure, plot(abs(dS_T2_analit),'g+--'), hold on,plot(abs(dS_T2_numer),'ro--')
% % abs(dS_T2_analit)./abs(dS_T2_numer);
% % figure, plot(abs(dS_T2_analit),'g+--'), hold on,plot(abs(dS_T2_numer)*0.0625,'ro--');
% % hold on
% % title('Analitical & Numberical for T2 derivatives')
% % legend('T2 analitical derivative','T2 numerical derivative')

%% 4 - Test EPG

% ... 4.1 - Parameters for generating signals (Obtain dictionary) ...
[All_dict_pp,All_pars_pp] = auxOptimiz_dict_pars_generator_ASM(T1_dic,T2_dic,B1_dic,...
    dTE,ETL,phase_refoc,phase_exc,FA_exc_dic, FA_refoc_dic,[],methodDic); % 10min

% ... 4.2 - get EPG with Shaihan toolbox
alpha0    = pi/180 * [90 repmat(flipAngle,ETL,1)'];  % vector of angles in radians
phase0    = pi/180 * [90 repmat(0,ETL,1)'];     % vector of phase in radians

% Get EPG
[F] = auxOptimiz_EPG_forward1(alpha0,phase0,'ESP',dTE,'T1',T1_dic,'T2',T2_dic); % set target, first simulate for ideal alpha
testF = F(1,3:end);

... 4.3 - Plots
% % % % figure,plot(abs(S),'ro--'), hold on,  plot(abs(All_dict_pp),'g+--'); hold on; plot(abs(testF),'bx--') 
% % % % hold on
% % % % title('Signal Mono-Exponential vs EPG')
% % % % legend('Signal Mono-Exponential','Signal my EPG','Signal EPG Shaihan')
% % % % 

fprintf('\n\n Sucessfully finished -  Dictionaries \n\n')


%% 5 - Get Gradient using Shaihan toolbox
% 5.1 - Parameters
Nch = 1;
Ny  = 2; Nx = Ny;    % TODO - resolution
B1  = ones(Nx*Ny,1); % TODO (add to function inputs) - Necessary shape Number of points - number of channels [Ns, Nch]

% 5.2 - Kronecker tensor product
frequencies = ones(1,ETL+1);                                                    % Frequency of each pulse
alpha_start = kron(ones(1,Nch),alpha0(cumsum(frequencies)));                    % Angle
phi_start   = kron(zeros(1,Nch),ones(1,length(alpha0(cumsum(frequencies)))));   % Phase
param_start = [real(alpha_start.*exp(1i*phi_start)), imag(alpha_start.*exp(1i*phi_start))];

initial        = 3;             % focus only on states from t=initial; 3 is first echo
Nt             = size(alpha_start,2) + 1;       % number of RF pulses (1x RF exc + ETL x RF refoc)
c              = zeros(1,Nt); 
c(initial:end) = 1; 

% 5.3 - get Objective function + Gradient
% % [obj_RF, grad_RF, FF_RF]  = auxOptimiz_obj_EPG13_RF(param_start,dTE,T1_dic,T2_dic,c,B1,frequencies,10e+06);
[obj_T2, grad_T2, FF_T2]  = auxOptimiz_obj_EPG13_T2(param_start,dTE,T1_dic,T2_dic,c,B1,frequencies,10e+06,valSigma);

% 5.4 - Get derivatives from myEPG
auxOptimiz_myEPGDerivativ

% 5.5 - Figures
% % figure()
% % subplot(121)
% % plot(abs(grad_RF(2:ETL+1)))
% % hold on
% % title('Gradient of EPG due to RF pulse')
% % subplot(122)
% % plot(abs(grad_T2(2:ETL+1)))
% % hold on
% % title('Gradient of EPG due to T2')

    % Plots of derivatives from mono-exponential expression
    
if testPlot == 'True'
    figure, plot(abs(dS_T2_analit),'ro--'), hold on,
    plot(abs(dS_T2_numer),'g+--'), hold on
    plot(abs(grad_T2(2:ETL+1)), 'bx--'), hold on
    title('T2 derivatives')
    legend('T2 analitical deriv.','T2 numerical deriv.','T2 EPG Shaihan deriv.')

    figure()
    subplot(151)
    plot(abs(dS_T2_analit),'ro--'), hold on,
    title('T2 analitical deriv.')
    subplot(152)
    plot(abs(dS_T2_numer),'g+--'), hold on
    title('T2 numerical deriv.')
    subplot(153)
    plot(abs(dS_dT2_myEPG),'c+--'), hold on
    title('T2 EPG num. deriv.')
    subplot(154)
    plot(abs(dS_dT2_Shahan),'y+--'), hold on
    title('T2 EPG num. Shaihan deriv.')
    subplot(155)
    plot(abs(grad_T2(2:ETL+1)), 'bx--'), hold on
    title('T2 EPG Shaihan ASM deriv.')

    aux_ratio = grad_T2(2:ETL+1)';
    ratio = dS_T2_analit ./ aux_ratio;
    figure();plot(ratio)
end

%% 6 - Comparison of uCRLB for mono-exponential vs EPG vs Shaihan

% === 6.1 - initial guess ============
% x0 = [dTE,ETL,flipAngle,]    
x0 = 16;

% === 6.2 - Set Objective ============
% [obj_T2, grad_T2, FF_T2]  = auxOptimiz_obj_EPG13_T2_test(param_start,dTE,T1_dic,T2_dic,c,B1,frequencies,10e+06,valSigma);
[obj_T2, grad_T2_Test, FF_T2]  = auxOptimiz_obj_EPG13_T2_test2(param_start,dTE,T1_dic,x0,c,B1,frequencies,10e+06,testPlot);

% === 6.3 - Plots of derivatives from mono-exponential expression =======
if testPlot == 'True'
    figure, plot(abs(dS_T2_analit),'ro--'), hold on,
    % plot(abs(dS_T2_numer),'g+--'), hold on
    % % plot(abs(grad_T2(2:ETL+1)), 'bx--'), hold on
    plot(abs(grad_T2_Test(2,ETL:end)),'g+--'), hold on
    title('T2 derivatives')
    % % legend('T2 analitical deriv.','T2 EPG Shaihan deriv.','T2 EPG Shaihan deriv. Test')
    legend('T2 analitical deriv.','T2 EPG Shaihan deriv. Test')

    % test gradient from derivative
    testGrad = abs(dS_T2_analit) ./ abs(grad_T2_Test(2,ETL:end));
    testGrad2 = abs(dS_T2_analit) ./ abs(dS_T2_numer);

    figure()
    plot(testGrad)
    hold on
    plot(testGrad2)
    legend('Relation between T2 analitical deriv. & dEPG','TRelation between T2 analitical deriv. & num deriv')

    % extra plots
    figure()
    subplot(151)
    plot(abs(dS_T2_analit),'ro--'), hold on,
    title('T2 analitical deriv.')
    subplot(152)
    plot(abs(dS_T2_numer),'g+--'), hold on
    title('T2 numerical deriv.')
    subplot(153)
    plot(abs(dS_dT2_myEPG),'c+--'), hold on
    title('T2 EPG num. deriv.')
    subplot(154)
    plot(abs(dS_dT2_Shahan),'y+--'), hold on
    title('T2 EPG num. Shaihan deriv.')
    subplot(155)
    plot(abs(grad_T2_Test(2,ETL:end))*(mean(testGrad)), 'bx--'), hold on
    title('T2 EPG labDerivativ.')
end

%%
% 6.3 - Building fmincon
constVariables.ETL = 30;
constVariables.dTE = 40;

% === 6.3.1 - Set Objective ============
% % [obj_T2, grad_T2, FF_T2]  = auxOptimiz_obj_EPG13_T2(param_start,dTE,T1_dic,T2_dic,c,B1,frequencies,10e+06,valSigma);

% === 6.3.2 - Define Parameters ========
    % 6.3.2.1 - initial guess
% x0 = [dTE,ETL,flipAngle,]    
x0 = 16;

% show initial objective
disp(['Initial Objective: ' num2str(obj_T2)])

    % 6.3.2.2 - variable bounds
% lb = [min(vector_dTE) min(vector_ETL) min(vector_flipAngle)]';  % lower-bounds
% ub = [max(vector_dTE) max(vector_ETL) max(vector_flipAngle)]';  % upper-bounds
lb = 1;  % lower-bounds
ub = 60;  % upper-bounds

    % 6.3.2.3 - linear constraints
A = [];
b = [];
Aeq = [];
beq = [];

    % 6.3.2.4 - nonlinear constraints
% % nonlincon = @(x) test_CRLB_con(x,constVariables);

    % 6.3.2.5 - options (exemplo Sbrizzi)
maxiter = 50;       % maximal iterations for algorithm
tolFun  = 1.0e-6;
MFE     = 125000;   % MaxFunEvals - exemplo FSE_optimization_ASM
options = optimoptions('fmincon','MaxIter',maxiter, 'MaxFunEvals',MFE, ...
    'GradObj','on','GradConstr','on','DerivativeCheck','off','TolFun',tolFun);  % options for the fmincon    
options = optimoptions(options,'Display','iter-detailed');
options = optimoptions(options,'PlotFcns',@optimplotfval);

% === 6.3.3 - optimize with fmincon ===
% % x = fmincon(objective,x0,A,b,Aeq,beq,lb,ub,nonlincon,options);
[x_fMinCon,fval] = fmincon(@(T2_test) auxOptimiz_obj_CRLB(param_start,dTE,T1_dic,T2_test,c,B1,frequencies,10e+06,valSigma),...
            x0,A,b,Aeq,beq,lb,ub,[],options);

% %     [sol1 fval] = fmincon(@(thet) obj_EPG13(thet,ESP,T1,T2,c,B1,TARGET,frequencies,klim),     ...
% %         param_start,[],[],[],[],[],[], ...
% %         @(thet) limit_peak_and_power(thet,max_power,max_a*pi,Nch,frequencies), ...
% %         options);
disp(['Optimized T2 value: ' num2str(x_fMinCon)])
disp(['Objective value is: ' num2str(fval)])

%% 


% % % % %% 6 - plots
% % % % sol = red2full(sol1,frequencies,Nch);
% % % % param_start1 = red2full(param_start,frequencies,Nch);
% % % % xopt = reshape(sol(1:Nt*Nch),Nt,Nch).'; %real parts
% % % % yopt = reshape(sol(Nt*Nch+1:end),Nt,Nch).'; % imaginary parts
% % % % xstart = reshape(param_start1(1:Nt*Nch),Nt,Nch).'; %real parts
% % % % ystart = reshape(param_start1(Nt*Nch+1:end),Nt,Nch).'; % imaginary parts
% % % % [obj, grad,FF1] = obj_EPG11(sol,ESP,T1,T2,c,B1,TARGET);
% % % % [obj, grad,FF] = obj_EPG11(param_start1,ESP,T1,T2,c,B1,TARGET);
% % % % N = size(FF,1)/2;
% % % % state = squeeze(FF1(2,:,:)+1i*FF1(N+2,:,:));
% % % % state0 = squeeze(FF(2,:,:)+1i*FF(N+2,:,:));
% % % % wTARGET = c*TARGET;
% % % % residual_abs = c*(abs(state)-abs(TARGET));
% % % % residual_compl = c*(state-TARGET);
% % % % residual_abs0 = c*(abs(state0)-abs(TARGET));
% % % % residual_compl0 = c*(state0-TARGET);
% % % % error_abs_optim = norm(residual_abs(:))/norm(wTARGET(:));
% % % % error_abs_start = norm(residual_abs0(:))/norm(wTARGET(:));
% % % % fprintf(1,'Absolute error before: %1.3f and after: %1.3f \n',...
% % % %     error_abs_start,error_abs_optim)
% % % % alpha_opt = abs(xopt+1i*yopt);
% % % % phi_opt = angle(xopt+1i*yopt);
% % % % Alpha_start = abs(xstart+1i*ystart);
% % % % Phi_start = angle(xstart+1i*ystart);
% % % % N = size(FF,1)/2;
% % % % 
% % % % 
% % % % %% 7 - Plots
% % % % figure;
% % % % subplot(2,1,1);stairs(180/pi*alpha_opt.','LineWidth',2);grid on;title('Tip angles, amplitude')
% % % % subplot(2,1,2);stairs(180/pi*phi_opt.','LineWidth',2);grid on;title('Tip angles, phases')
% % % % figure;
% % % % subplot(5,2,2);
% % % % plot(180/pi*alpha_opt.','LineWidth',2);ylim([0 180*max_a]);
% % % % title('OPTIMIZED \alpha');
% % % % subplot(5,2,4);
% % % % plot(180/pi*phi_opt.','LineWidth',2);ylim([-180 180]);
% % % % title('OPTIMIZED \phi');
% % % % subplot(5,2,1);
% % % % plot(180/pi*Alpha_start.','LineWidth',2);ylim([0 180*max_a]);
% % % % title('STARTING \alpha');
% % % % subplot(5,2,3);
% % % % plot(180/pi*Phi_start.','LineWidth',2);ylim([-180 180]);
% % % % title('STARTING \phi');
% % % % subplot(5,2,5);
% % % % plot(squeeze(abs((FF(2,:,:)+1i*FF(N+2,:,:)))),'LineWidth',.5);title('STARTING SIGNAL');xlim([1 Nt+1]);grid on;
% % % % hold on;plot(abs(sum(TARGET,2)/Ns),'bo','LineWidth',2);
% % % % subplot(5,2,6);
% % % % plot(squeeze(abs((FF1(2,:,:)+1i*FF1(N+2,:,:)))),'LineWidth',.5);title('OPTIMIZED SIGNAL');xlim([1 Nt+1]);grid on;
% % % % hold on;plot(abs(sum(TARGET,2)/Ns),'bo','LineWidth',2);
% % % % subplot(5,2,7);
% % % % imagesc(abs(squeeze(FF(2,:,:)+1i*FF(N+2,:,:))).');title('STARTING SIGNAL');grid on;ylabel('voxels')
% % % % subplot(5,2,8);
% % % % imagesc(abs(squeeze(FF1(2,:,:)+1i*FF1(N+2,:,:))).');title('OPTIMIZED SIGNAL');grid on;ylabel('voxels')
% % % % subplot(5,2,9);
% % % % imagesc(angle(squeeze(FF(2,:,:)+1i*FF(N+2,:,:))).');title('STARTING SIGNAL');grid on;ylabel('voxels')
% % % % subplot(5,2,10);
% % % % imagesc(angle(squeeze(FF1(2,:,:)+1i*FF1(N+2,:,:))).');title('OPTIMIZED SIGNAL');grid on;ylabel('voxels')

a = 1

