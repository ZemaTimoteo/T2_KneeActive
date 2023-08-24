% Main script for optimal control EPG design - Sbrizzi paper - TEST
%
% by Alessandro Sbrizzi, June 2015 (email: a.sbrizzi@umcutrecht.nl)
% based on EPG simulator implementation of Shaihan Malik
% adapted TTFernandes, 2023

%% 0 - Get Path's
clear all
close all
cd('D:\Tiago\Trabalho\2021_2025_PhD\Projects\qMRI_Joint\Code\matlabCode\Toolboxes\optimal-control-EPG-master')

%% 1 - read info data (two examples are given here)
info_test4;
%info_test5;

%% 2 - Parameters
% ... 2.1 - load B1 map ...
% % load(B1_file); % pode ser = 1 desde que corresponde a dimensoes
B1          = 1;
Nch         = size(B1,2);
frequencies = [ones(1,singles), const*ones(1,intervals)];
figure;bar(frequencies);title('tip angles repetitions (frequencies)');grid on;
Nt   = sum(frequencies);% number of time points
klim = min([klim, round(Nt/2)]);


%% 3 - Target 
% ... 3.1 - Set angles for target ...
alpha0         = transpose(a0(:));
c              = zeros(1,Nt+1);
c(initial:end) = 1; 

% ... 3.2 - Set target, first simulate for ideal alpha ...
[F]    = EPG_forward1(alpha0,'ESP',ESP,'T1',T1,'T2',T2); 
target = full(F(2,:)); % target F vector (containes Nt-1 echoes)
target = abs(target);
Ns     = size(B1,1);

% ... 3.3 - Starting value. Start with quadrature values ...
alpha_start = kron(ones(1,Nch),alpha0(cumsum(frequencies)));                    % Kronecker tensor product
phi_start   = kron(zeros(1,Nch),ones(1,length(alpha0(cumsum(frequencies)))));   % Kronecker tensor product
TARGET      = repmat(target(:),[ 1 Ns]);                                        % target for all spatial locations
param_start = [real(alpha_start.*exp(1i*phi_start)), imag(alpha_start.*exp(1i*phi_start))];

%% 4 - Optimization

% ... 4.1 - Parameters ...
param_all      = red2full(param_start,frequencies,Nch);
max_power      = pow_factor*norm(param_all(:))^2;                                 % constraint power to power of starting sequence
max_chpow      = chpow_factor*ones(1,Nch)*sum(alpha0.^2);

% ... 4.2 - Cost Function ...
[obj, grad,FF] = obj_EPG13(param_start,ESP,T1,T2,c,B1,TARGET,frequencies,10e+06); % funcao de custo que quero minimizar
N       = size(FF,1)/2;

% ... 4.3 - Options ...
options = optimoptions('fmincon','MaxIter',maxiter, 'MaxFunEvals',1000*length(alpha_start)*5, ...
    'GradObj','on','GradConstr','on','DerivativeCheck','off','TolFun',1.0e-09); % minimo da f. de custo mas com constrains
options = optimoptions(options,'Display','iter-detailed'); % plot and display each iteration:
options = optimoptions(options,'PlotFcns',@optimplotfval); % plot and display each iteration:

% ... 4.4 - Optimization ...
tic

%%% Phase average over time
Z_TARGET           = squeeze(angle(FF(2,:,:)+1i*FF(N+2,:,:)));
Z_TARGETs          = unwrap(Z_TARGET,[],1);
Z_TARGETs(2:end,:) = repmat(mean(Z_TARGETs(initial:end,:),1),[Nt 1]);
TARGET             = abs(TARGET).*exp(1i*Z_TARGETs);

switch pow_constr % constrains
    case 0
        %       fmincon(FUN                                                           ,X          ,A ,B ,Aeq,Beq,LB,UB,NONLCON                          ,options)
        
        
        % x =   fmincon(objective                                                     ,x0         ,A ,b ,Aeq,beq,lb,ub,nonlincon                        ,options);
        %       fmincon(@(thet) obj_EPG13(thet,ESP,T1,T2,c,B1,TARGET,frequencies,klim),param_start,[],[],[] ,[] ,[],[],@(thet) limit_RF(thet,1.5*pi,Nch),options);         
    % matlab implementation:
    %[sol1 fval] = fmincon(@(thet) obj_EPG13(thet,ESP,T1,T2,c,B1,TARGET,frequencies,klim),param_start,[],[],[],[],[],[],@(thet) limit_RF(thet,1.5*pi,Nch),options); 
    % c++ implementation:
    [sol1 fval] = fmincon(@(thet) obj_EPG13_cpp(thet,ESP,T1,T2,c,B1,TARGET,frequencies,klim),param_start,[],[],[],[],[],[],@(thet) limit_RF(thet,max_a*pi,Nch),options); 
    case 1
    % c++ implementation with constr on peak AND total power
    [sol1 fval] = fmincon(@(thet) obj_EPG13_cpp(thet,ESP,T1,T2,c,B1,TARGET,frequencies,klim), ...
    param_start,[],[],[],[],[],[],@(thet) limit_peak_and_power(thet,max_power,max_a*pi,Nch,frequencies),options); 
    case 2
    % c++ implementation with constr on peak AND total CHANNEL power
    [sol1 fval] = fmincon(@(thet) obj_EPG13_cpp(thet,ESP,T1,T2,c,B1,TARGET,frequencies,klim), ...
    param_start,[],[],[],[],[],[],@(thet) limit_peak_and_ch_power(thet,max_chpow,max_a*pi,Nch,frequencies),options); 
end
N = size(FF,1)/2;
toc

%% 5 - Plots1
sol = red2full(sol1,frequencies,Nch);
param_start1 = red2full(param_start,frequencies,Nch);
xopt = reshape(sol(1:Nt*Nch),Nt,Nch).'; %real parts
yopt = reshape(sol(Nt*Nch+1:end),Nt,Nch).'; % imaginary parts
xstart = reshape(param_start1(1:Nt*Nch),Nt,Nch).'; %real parts
ystart = reshape(param_start1(Nt*Nch+1:end),Nt,Nch).'; % imaginary parts
[obj, grad,FF1] = obj_EPG11(sol,ESP,T1,T2,c,B1,TARGET);
[obj, grad,FF] = obj_EPG11(param_start1,ESP,T1,T2,c,B1,TARGET);
N = size(FF,1)/2;
state = squeeze(FF1(2,:,:)+1i*FF1(N+2,:,:));
state0 = squeeze(FF(2,:,:)+1i*FF(N+2,:,:));
wTARGET = c*TARGET;
residual_abs = c*(abs(state)-abs(TARGET));
residual_compl = c*(state-TARGET);
residual_abs0 = c*(abs(state0)-abs(TARGET));
residual_compl0 = c*(state0-TARGET);
error_abs_optim = norm(residual_abs(:))/norm(wTARGET(:));
error_abs_start = norm(residual_abs0(:))/norm(wTARGET(:));
fprintf(1,'Absolute error before: %1.3f and after: %1.3f \n',...
    error_abs_start,error_abs_optim)
alpha_opt = abs(xopt+1i*yopt);
phi_opt = angle(xopt+1i*yopt);
Alpha_start = abs(xstart+1i*ystart);
Phi_start = angle(xstart+1i*ystart);
N = size(FF,1)/2;

%% 6 - Plots2
figure;
subplot(2,1,1);stairs(180/pi*alpha_opt.','LineWidth',2);grid on;title('Tip angles, amplitude')
subplot(2,1,2);stairs(180/pi*phi_opt.','LineWidth',2);grid on;title('Tip angles, phases')
figure;
subplot(5,2,2);
plot(180/pi*alpha_opt.','LineWidth',2);ylim([0 180*max_a]);
title('OPTIMIZED \alpha');
subplot(5,2,4);
plot(180/pi*phi_opt.','LineWidth',2);ylim([-180 180]);
title('OPTIMIZED \phi');
subplot(5,2,1);
plot(180/pi*Alpha_start.','LineWidth',2);ylim([0 180*max_a]);
title('STARTING \alpha');
subplot(5,2,3);
plot(180/pi*Phi_start.','LineWidth',2);ylim([-180 180]);
title('STARTING \phi');
subplot(5,2,5);
plot(squeeze(abs((FF(2,:,:)+1i*FF(N+2,:,:)))),'LineWidth',.5);title('STARTING SIGNAL');xlim([1 Nt+1]);grid on;
hold on;plot(abs(sum(TARGET,2)/Ns),'bo','LineWidth',2);
subplot(5,2,6);
plot(squeeze(abs((FF1(2,:,:)+1i*FF1(N+2,:,:)))),'LineWidth',.5);title('OPTIMIZED SIGNAL');xlim([1 Nt+1]);grid on;
hold on;plot(abs(sum(TARGET,2)/Ns),'bo','LineWidth',2);
subplot(5,2,7);
imagesc(abs(squeeze(FF(2,:,:)+1i*FF(N+2,:,:))).');title('STARTING SIGNAL');grid on;ylabel('voxels')
subplot(5,2,8);
imagesc(abs(squeeze(FF1(2,:,:)+1i*FF1(N+2,:,:))).');title('OPTIMIZED SIGNAL');grid on;ylabel('voxels')
subplot(5,2,9);
imagesc(angle(squeeze(FF(2,:,:)+1i*FF(N+2,:,:))).');title('STARTING SIGNAL');grid on;ylabel('voxels')
subplot(5,2,10);
imagesc(angle(squeeze(FF1(2,:,:)+1i*FF1(N+2,:,:))).');title('OPTIMIZED SIGNAL');grid on;ylabel('voxels')
