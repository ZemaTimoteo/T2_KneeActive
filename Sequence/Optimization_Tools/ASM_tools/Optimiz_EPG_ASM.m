
% ... 4.1 - Parameters for generating signals with slr_profile.py and my_epg.py ...
% Obtain dictionary
[All_dict_pp{ii},All_pars_pp{ii}] = auxOptimiz_dict_pars_generator_ASM(T1_dic,T2_dic,B1_dic,...
    dTE,ETL,phase_refoc,phase_exc,FA_exc_dic, FA_refoc_dic,file_path_data,methodDic); % 10min
fprintf(['      Successfull models with obtained dictionary ',num2str(ii),' / ',...
    num2str(size(Models_accepted,1)),'\n'])

%Models_accepted: [TE(ms) #ETL angle Nx/Ny TimeScan(min) b1+rms(uT) Trec(ms) TR_scan(s)]
fprintf('\n\n Sucessfully finished -  Dictionaries \n\n')

