%% Parameters
gg = 1;

%% Signal 8
ii = 8;
All_dict = All_dict_pp{ii};
All_pars = All_pars_pp{ii};
% Signals
S_8  = All_dict(:,find(All_pars(:,1) == vector_T2(gg) & All_pars(:,3) == B1)); % Signal S - Dictionary
% TRvar
auxVariab = exp( - Models_accepted(ii,7)/T1_dic);
aux2var   = (1 - auxVariab)/ (1 - auxVariab*S_8(Models_accepted(ii,2)))
S_8       = aux2var.* S_8;

%% Signal 60
% ii = 60;
% All_dict = All_dict_pp{ii};
% All_pars = All_pars_pp{ii};
% % Signals
% S_60  = All_dict(:,find(All_pars(:,1) == vector_T2(gg) & All_pars(:,3) == B1)); % Signal S - Dictionary
% % TRvar
% auxVariab = exp( - Models_accepted(ii,7)/T1_dic);
% aux2var   = (1 - auxVariab)/ (1 - auxVariab*S_60(Models_accepted(ii,2)))
% S_60      = aux2var.* S_60;

%% Signal 156
ii = 96;
All_dict = All_dict_pp{ii};
All_pars = All_pars_pp{ii};
% Signals
S_156  = All_dict(:,find(All_pars(:,1) == vector_T2(gg) & All_pars(:,3) == B1)); % Signal S - Dictionary
% TRvar
auxVariab = exp( - Models_accepted(ii,7)/T1_dic);
aux2var   = (1 - auxVariab)/ (1 - auxVariab*S_156(Models_accepted(ii,2)))
S_156     = aux2var.* S_156;

%% Plots
figure()
subplot(121)
plot(abs(S_8),'r')
hold on
% plot(abs(S_60),'g')
% hold on
% legend('ETL60')
plot(abs(S_156),'b')
hold on
legend('ETL8','ETL60','ETL156')
subplot(122)
plot(abs(S_8-S_156(1:6)))
% % hold on
% % plot(abs(S_8-S_60(1:6)))
