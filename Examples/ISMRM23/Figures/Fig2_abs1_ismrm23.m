%% ISMRM23 Figure2 abs1
clear all
cd('Github\Examples\ISMRM23')
load('originalValues_tests.mat')

%% Get data
% T2 values
TV     = Table1S1(:,1); % True values
T2_200 = Table1S1(:,2); 
T2_45  = Table1S1(:,4);
T2_8   = Table1S1(:,6);

% error associated
errT2_200 = Table1S1(:,3);
errT2_45  = Table1S1(:,5);
errT2_8   = Table1S1(:,7);

% (%) of how close to real value
percT2_200 = (abs(TV-T2_200))./TV*100;
percT2_45 = (abs(TV-T2_45))./TV*100;
percT2_8 = (abs(TV-T2_8))./TV*100;

% error associated
err_percT2_200 = Table1S1(:,3)./TV*100;
err_percT2_45  = Table1S1(:,5)./TV*100;
err_percT2_8   = Table1S1(:,7)./TV*100;

%% Value closeness
% % T2est = figure()
% % plot(TV,'r','linewidth',2)
% % hold on
% % errorbar(T2_200,errT2_200)
% % errorbar(T2_45,errT2_45)
% % errorbar(T2_8,errT2_8)
% % hold on
% % title('True value VS Estimations','FontSize',20)
% % xlabel('Vials','FontSize',15)
% % ylabel('T2 (ms)','FontSize',15)
% % ylim([-5 650])
% % xlim([0.4 14.5])
% % % % plot([-0.5 14.5],[0 0],'k:','LineWidth',1)
% % legend('True Values','T_2 opt. 200ms','T_2 opt. 45ms','T_2 opt. 8ms','Location','northeast','FontSize',15)
% % xticks([0:14])
% % 
% % saveas(T2est,['Fig2a_abs1_ismrm23.png'])

%% Percentage plot
% % percent_erro = figure()
% % % % errorbar(percT2_200,log(err_percT2_200),'linewidth',2)
% % errorbar(percT2_200,log(err_percT2_200),'.','linewidth',2,'MarkerSize',12)
% % hold on
% % % % errorbar(percT2_45,log(err_percT2_45),'linewidth',2)
% % % % errorbar(percT2_8,log(err_percT2_8),'linewidth',2)
% % errorbar(percT2_45,log(err_percT2_45),'.','linewidth',2,'MarkerSize',12)
% % errorbar(percT2_8,log(err_percT2_8),'.','linewidth',2,'MarkerSize',12)
% % title('% Bias relative to Reference','FontSize',20)
% % xlabel('Vials','FontSize',15)
% % ylabel('% of True value','FontSize',15)
% % ylim([-5 105])
% % xlim([0.4 14.5])
% % plot([-0.5 14.5],[0 0],'k:','LineWidth',1.5)
% % plot([-0.5 14.5],[10 10],'r:','LineWidth',1.5)
% % legend('T_2 opt. 200ms','T_2 opt. 45ms','T_2 opt. 8ms','Zero Line','10% Line','Location','northwest','FontSize',15)
% % xticks([0:14])
% % 
% % saveas(percent_erro,['Fig2b_abs1_ismrm23.png'])

%% FIg 2 a - test - Value closeness

T2est = figure()

ylowerBound = -5;
yUpperBound = 750;

x0  = 1:7:size(percT2_200,1)*7;
x1  = x0+1;
x2  = x0+2;
x3  = x0+3;
x4  = x0+4;
x5  = x0+5;

scatter(x1,TV,'r','linewidth',1.5)
hold on
errorbar(x4,T2_200,(errT2_200),'.','MarkerSize',12,'linewidth',1.5)
errorbar(x3,T2_45,(errT2_45),'.','MarkerSize',12,'linewidth',1.5)
errorbar(x2,T2_8,(errT2_8),'.','MarkerSize',12,'linewidth',1.5)

hold on
title('True value VS Estimations','FontSize',20)
xlabel('Vials','FontSize',15)
ylabel('T2 (ms)','FontSize',15)
ylim([ylowerBound yUpperBound])
xlim([0.4 97.5])
% % plot([-0.5 14.5],[0 0],'k:','LineWidth',1)
xticks([3:7:97])
xticklabels({'1','2','3','4','5','6','7','8','9','10','11','12','13','14'})
for i=1:14
    xx = [x0(i) x5(i) x5(i) x0(i)];
    yy = [ylowerBound ylowerBound yUpperBound yUpperBound];
    patch(xx,yy,'green','FaceAlpha',.03)
end
legend('True Values','T_2 opt. 200ms','T_2 opt. 45ms','T_2 opt. 8ms','Location','northeast','FontSize',15)


saveas(T2est,['Fig2a_NEW_abs1_ismrm23.png'])


%% FIg 2 b - test - Percentage plot

ylowerBound = -1;
yUpperBound = 105;

x0  = 1:6:size(percT2_200,1)*6;
x1  = x0+1;
x2  = x0+2;
x3  = x0+3;
x4  = x0+4;

percent_erro = figure()
errorbar(x3,percT2_200,(err_percT2_200),'.','MarkerSize',12,'linewidth',1.5)
hold on
errorbar(x2,percT2_45,(err_percT2_45),'.','MarkerSize',12,'linewidth',1.5)
errorbar(x1,percT2_8,(err_percT2_8),'.','MarkerSize',12,'linewidth',1.5)

% plot(x4,zeros1,'b','.')
ylim([ylowerBound yUpperBound])

hold on
for i=1:14
    xx = [x0(i) x4(i) x4(i) x0(i)];
    yy = [ylowerBound ylowerBound yUpperBound yUpperBound];
    patch(xx,yy,'blue','FaceAlpha',.05)
end

title('% Bias relative to Reference','FontSize',20)
xlabel('Vials','FontSize',15)
ylabel('% of True value','FontSize',15)
ylim([-5 105])
plot([-0.5 83.5],[0 0],'k:','LineWidth',1.5)
plot([-0.5 83.5],[10 10],'r:','LineWidth',1.5)
legend('T_2 opt. 200ms','T_2 opt. 45ms','T_2 opt. 8ms','Location','northwest','FontSize',15)
xticks([3:6:83])
xticklabels({'1','2','3','4','5','6','7','8','9','10','11','12','13','14'})

saveas(percent_erro,['Fig2b_NEW_abs1_ismrm23.png'])
