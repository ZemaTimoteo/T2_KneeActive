% Read T2 maps for Repeatability study from Siemens Dicom files and export
% maps
%
% @TTFernandes, November, 2022

%% 0 - Set matlab paths and toolboxes
clear all
clc

addpath(genpath('D:\Tiago\Trabalho\2021_2025_PhD\Projects\qMRI_Joint\Code\matlabCode\MSE_recon\aux_Functions'))
addpath(genpath('D:\Tiago\Trabalho\2021_2025_PhD\Tutoria\2021-2022_3P_Estagio\Code\matlab\Dictionary_aux_Functions'));
addpath(genpath('D:\Tiago\Trabalho\2021_2025_PhD\Projects\GIRFOS_tool\Code\0_matlabCode\10_toolboxes'))
addpath(genpath('D:\Tiago\Trabalho\2021_2025_PhD\Projects\qMRI_Joint\Code\matlabCode\MSE_recon'))

%% 0.5 - Tests
plotTest     = 'True';     % 'Fals' or 'True'
saveResults  = 'True';     % 'Fals' or 'True'
maskTest     = 'ready';    % 'getIt' - Get mask OR 'ready' - load mask
avgBef_T2est = 'Fals';     % 'True' if avg of data before T2 estimation or 'Fals'

sliceTest = [1];
numSli    = sliceTest;   % Number of Slices

%% 1 - LOAD Data% select directory for recon.
myCD = 'C:\Users\filia\Documents\PhD\Projetos\qMRI\Sequences\MSE\RepeatabilityStudy_ISMRM23\vMSE';
cd(myCD)

% get dir / select data from day and run folder
str = {'select file'};
data_dir=uigetdir;
cd(data_dir)
data_listing = dir(data_dir);

% --  Get Dir_results
dir_results = [data_dir '\results'];
mkdir(dir_results)

% load file
for i=3:size(data_listing,1)-1
    im_data_file   = [data_dir filesep data_listing(i).name];   % Dicom Data
    recon_Img_full = (dicomread(im_data_file));
    infoDicom9     = (dicominfo(im_data_file));
    image{i-2}     = recon_Img_full;
end

% Get name
vT2dic_name = ['vT2_map_ data_' data_dir(end-10:end-7)];

% Ger imSE_final 
imgSE_final = zeros(size(image{1},1),size(image{1},2),size(image,2));
for j=1:size(image,2)
    imgSE_final(:,:,j)=double(image{j});
end

clear data_listing image recon_Img_full

%% 1 - load Data Phantom MSE
% ... 1.1 - select directory for recon ...
cd(data_dir)

% ... 1.2 - get parameters ...
Ny        = size(imgSE_final,2);
nTE       = infoDicom9.EchoNumbers;     % Number of Echos of TEs
ESP       = infoDicom9.EchoTime*1e-4;   % First TE - in (s)
TE_first  = ESP;                        % First TE - in (s)
TE_vector = [ESP:ESP:nTE*ESP];          % TE vector - in (s)
nl        = Ny;                         % Image size
slice_inf = min(sliceTest);
slice_sup = max(sliceTest);
t         = TE_vector';

% figures
figure()
for ii=1:size(imgSE_final,3)
    subplot(2,5,ii)
    imshow(imgSE_final(:,:,ii)./max(max(max(imgSE_final))),[])     
end
fprintf('\n\n 1 - Sucessfully finished - Load\n\n')


%% 2 - inputs

if maskTest == 'ready'
    nvials = 14;
    for jj=1:nvials
        % --  Circular Vial mask
        recon_Img_full = imgSE_final(:,:,1);
        T2_mask_m = zeros(size(recon_Img_full));
        figure('Name','Results','WindowState','maximized'); imshow(recon_Img_full,[]);           % display image 1st echo random slice with good cartilage display
        title(['Create Mask for Preprocessing vial=', num2str(jj)])
        
        circle = drawcircle('Color','r');  % select roi with cartilage region only
        r      = round(circle.Radius);     % circle radius
        c_x0   = round(circle.Center(1));  % x0 position
        c_y0   = round(circle.Center(2));  % y0 position
        
        disc            = filldisc(r);       % fill circle
        [disc_x,disc_y] = size(disc);
        close
        T2_mask_m(c_y0-r:c_y0+r,c_x0-r:c_x0+r) = disc;
        T2_mask(:,:,jj) = double(T2_mask_m);
        
        % -- Save masks
        cd(dir_results)        
        name_mask = 'roiMSE_T2Phantom_py.mat';
        save(name_mask,'T2_mask')       
    end
    
else    % -- Load masks
    cd(dir_results)
    name_mask = 'roiMSE_T2Phantom_py.mat';  
    load(name_mask)    
end

fprintf('\n\n 2 - Sucessfully finished - Get Inputs\n\n')


%% 3 - T2 mapping from mono-exponential fit (excluding 1st echo)
if avgBef_T2est == 'Fals'
    tic

    if saveResults == 'True'
        T2_exp_map            = zeros(nl,nl,numSli);
        T2_wt1echo_exp_map    = zeros(nl,nl,numSli);
        T2_evenecho_exp_map   = zeros(nl,nl,numSli);
        ShortT2_biexp_wt1_map = zeros(nl,nl,numSli);
        LongT2_biexp_wt1_map  = zeros(nl,nl,numSli);

        n_vial = 1;

        for ind = 1:nvials 
            X_SE = squeeze(imgSE_final(:,:,:));
            X_SE = permute(X_SE,[3 2 1]);
            X_SE = reshape(X_SE,[nTE nl*nl]);
            X(:,:) = X_SE;  % Data image matrix
            aux_T2_mask    = reshape(T2_mask(:,:,ind),[1 size(T2_mask,1)*size(T2_mask,2)]);

            for i = 1:size(X,2)
                if aux_T2_mask(1,i) == 0

                    [T2_exp(1,i)]         = 0;
                    [T2_wt1echo_exp(1,i)] = 0;
                    [T2_even_exp(1,i)]    = 0;
                    [shortT2(1,i)]        = 0;
                    [longT2(1,i)]         = 0;

                    fit_out{ind,i}         = 0;
                    fit_wt1echo_out{ind,i} = 0;
                    fit_even_out{ind,i}    = 0;
                else

                    ydat = squeeze(abs(X(:,i)));
                    if sum(ydat) == 0
                        [T2_exp(1,i)]         = 0;
                        fit_out{ind,i}        = 0;
                    else               
                        [T2_exp(1,i) , fit_out{ind,i}]                 = monoExpfit(t,double(ydat));
                    end
                    if nTE>2
                        if sum(double(ydat(2:end))) == 0                
                            [T2_wt1echo_exp(1,i)]  = 0;
                            fit_wt1echo_out{ind,i} = 0;
                        else        
                            [T2_wt1echo_exp(1,i) , fit_wt1echo_out{ind,i}] = monoExpfit(t(2:end),double(ydat(2:end)));                    
                        end
                        if sum(double(ydat(2:2:end))) == 0
                            [T2_even_exp(1,i)]  = 0;                    
                            fit_even_out{ind,i} = 0;
                        else
                            [T2_even_exp(1,i) , fit_even_out{ind,i}]       = monoExpfit(t(2:2:end),double(ydat(2:2:end)));
                        end
                    end
                end
                fprintf(['   ----  Mono-exponential fit for Point ', num2str(i),' / ',num2str(size(X,2)),'   ------ \n'])

            end    

            T2_exp_map (:,:,n_vial) = reshape(abs(T2_exp),[nl,nl]);

            if nTE>2        
                T2_wt1echo_exp_map (:,:,n_vial)  = reshape(abs(T2_wt1echo_exp),[nl,nl]);        
                T2_evenecho_exp_map (:,:,n_vial) = reshape(abs(T2_even_exp),[nl,nl]);
            end

            disp(strcat('T2 mono exp fit N slice = ',num2str(n_vial)));

            n_vial = n_vial + 1;

        end

        % --- save results ---
        cd(data_dir)    
        if nTE>2
            save(['resultsT2_monoexp_maps.mat'],'T2_exp_map','T2_wt1echo_exp_map','T2_evenecho_exp_map')
        else
            save(['resultsT2_monoexp_maps.mat'],'T2_exp_map')
        end
    else
        load('resultsT2_monoexp_maps.mat')
    end

    % ...means ...
    mean_T2_exp_map          = mean(mean(nonzeros(T2_exp_map)));          % in ms
    if nTE>2
        mean_T2_wt1echo_exp_map  = mean(mean(nonzeros(T2_wt1echo_exp_map)));  % in ms
        mean_T2_evenecho_exp_map = mean(mean(nonzeros(T2_evenecho_exp_map))); % in ms
    end

    % .. get maxValue ..
    maxVal = max(max(max(max((double(abs(T2_exp_map)))))));
    maxVal = 1000;

    % .. Plots ..
    if nTE<=2
        figure()
        montage_T2map_exp = permute(T2_exp_map,[1 2 4 3]);
        montage(montage_T2map_exp,[0 maxVal]);hold on; colormap hot;
        title(['T2 map mono exponential = ' num2str(mean_T2_exp_map)]);
        fprintf(['T2 map mono exponential = ' num2str(mean_T2_exp_map) 'ms \n']);

    else
        figure()
        montage_T2map_exp = permute(T2_exp_map,[1 2 4 3]);
% %         subplot(131); 
        montage(montage_T2map_exp,[0 maxVal]);hold on; colormap hot;
        title(['T2 map mono exponential = ' num2str(mean_T2_exp_map)]);
        fprintf(['T2 map mono exponential = ' num2str(mean_T2_exp_map) 'ms \n']);

        montage_T2map_exp_wt1echo = permute(T2_wt1echo_exp_map,[1 2 4 3]);
        subplot(132);montage(montage_T2map_exp_wt1echo,[0 maxVal]);colormap hot;
        title(['T2 map mono exponential excluding 1st echo = ' num2str(mean_T2_wt1echo_exp_map)]);
        fprintf(['T2 map mono exponential excluding 1st echo = ' num2str(mean_T2_wt1echo_exp_map) 'ms \n']);

        montage_T2map_exp_evenecho = permute(T2_evenecho_exp_map,[1 2 4 3]);
        subplot(133); montage(montage_T2map_exp_evenecho,[0 maxVal]);colormap hot;
        title(['T2 map mono exponential even echoes = ' num2str(mean_T2_evenecho_exp_map)]);
        fprintf(['T2 map mono exponential even echoes = ' num2str(mean_T2_evenecho_exp_map) 'ms \n']);
    end

    toc


    fprintf('\n\n 3 - Sucessfully finished - T2 mapping from mono-exponential fit\n\n')

    %% test vials T2 values
    nvials    = 1; % number of vials
    Vial_mask = zeros(size(montage_T2map_exp,1),size(montage_T2map_exp,1),nvials);

    %manually entries
    % --- B) - 128by128 ----------------------------------------------------------------
    v1_iniY = 72; v1_endY = 75; v1_iniX = 74; v1_endX = 78; 
    % ----------------------------------------------------------------------------------

    Vial_mask(v1_iniX:v1_endX-1,v1_iniY:v1_endY-1,1) = ones(size(1:v1_endX-v1_iniX,2),size(1:v1_endY-v1_iniY,2));

    % % figure(); imshow(Vial_mask(:,:,2),[0 max(max(Vial_mask(:,:,2)))]); %display image 1st echo random slice with good cartilage display
    % % figure(); imshow(Vial_mask(:,:,1),[0 max(max(Vial_mask(:,:,1)))]); %display image 1st echo random slice with good cartilage display

    for nv = 1:nvials
        aux_meanVial_T2map_exp(:,:,nv) = Vial_mask(:,:,nv).*montage_T2map_exp;
        meanVial_T2map_exp(nv)         = mean(nonzeros( aux_meanVial_T2map_exp(:,:,nv)));

        if nTE>2
            aux_meanVial_T2map_exp_wt1echo(:,:,nv)  = Vial_mask(:,:,nv).*montage_T2map_exp_wt1echo;
            aux_meanVial_T2map_exp_evenecho(:,:,nv) = Vial_mask(:,:,nv).*montage_T2map_exp_evenecho;

            meanVial_T2map_exp_wt1echo(nv)  = mean(nonzeros( aux_meanVial_T2map_exp_wt1echo(:,:,nv)));
            meanVial_T2map_exp_evenecho(nv) = mean(nonzeros( aux_meanVial_T2map_exp_evenecho(:,:,nv)));      
        end
    end

    % print values
    T2map   = {'Theoric';'T2map_exp';'T2map_exp_wt1echo';'T2map_exp_evenecho'};
    vial_1  = [315;meanVial_T2map_exp(1);meanVial_T2map_exp_wt1echo(1);meanVial_T2map_exp_evenecho(1)];
    T = table(T2map,vial_1)
    
    
    %% for test all vials
    clear Vial_mask
    nvials    = 1;
    Vial_mask = zeros(size(montage_T2map_exp,1),size(montage_T2map_exp,1),nvials);
    
    % --- B) - get poinst for vials 128by128
    v1_iniY  = 37;  v1_endY  = 40;      v1_iniX  = 63; v1_endX  = 67;
    
    % get masks for vials
    Vial_mask(v1_iniX:v1_endX-1,v1_iniY:v1_endY-1,1)      = ones(size(1:v1_endX-v1_iniX,2),size(1:v1_endY-v1_iniY,2));

% %     figure(); imshow(Vial_mask(:,:,9),[0 max(max(Vial_mask(:,:,9)))]); %display image 1st echo random slice with good cartilage display

    for nv = 1:nvials
        aux_meanVial_T2map_exp(:,:,nv) = Vial_mask(:,:,nv).*montage_T2map_exp;
        meanVial_T2map_exp(nv)         = mean(nonzeros( aux_meanVial_T2map_exp(:,:,nv)));
        stdVial_T2map_exp(nv)          = std(nonzeros( aux_meanVial_T2map_exp(:,:,nv)));

        if nTE>2
            aux_meanVial_T2map_exp_wt1echo(:,:,nv)  = Vial_mask(:,:,nv).*montage_T2map_exp_wt1echo;
            aux_meanVial_T2map_exp_evenecho(:,:,nv) = Vial_mask(:,:,nv).*montage_T2map_exp_evenecho;

            meanVial_T2map_exp_wt1echo(nv)  = mean(nonzeros( aux_meanVial_T2map_exp_wt1echo(:,:,nv)));
            meanVial_T2map_exp_evenecho(nv) = mean(nonzeros( aux_meanVial_T2map_exp_evenecho(:,:,nv)));      
        end
    end

    % print values
    T2map       = {'Theoric';'T2map_exp';'T2map_exp_wt1echo';'T2map_exp_evenecho'};
    vial_1      = [581.33;meanVial_T2map_exp(1);meanVial_T2map_exp_wt1echo(1);meanVial_T2map_exp_evenecho(1)];

    T = table(T2map,vial_1)

    % save Table            
    save(['table_of_T2values.mat'],'T')

end

%% 5 - Average before T2 maps
if avgBef_T2est == 'True'
    
    % 5.0 - clear variables
    clear T2_exp_map fit_out meanVial_Intensity aux_meanVial
    nvials    = 1;
    
    % 5.1 - Get all vials
    clear Vial_mask
    Vial_mask  = zeros(size(imgSE_final,1),size(imgSE_final,1),nvials);
    Vial_noise = zeros(size(imgSE_final,1),size(imgSE_final,1));
    SNR        = zeros(nTE,nvials);
    theoric_T2Vial = [315];    
    
    % 5.1.5 - Echo_Start the analysis
    echoStart    = 2;
    echoEnd      = zeros(1,nvials);
    threh_TEvect = theoric_T2Vial*1e-3*5;
    
    for nv = 1:nvials        
        aux1 = (threh_TEvect(nv)>=TE_vector);
        try
            aux_echoEnd = find( aux1 == 0 );
            echoEnd(nv) = aux_echoEnd(1);
        catch
            echoEnd(nv) = size(imgSE_final,4);
        end
    end
            
% %     vNoise_iniY = 11; vNoise_endY = 55;     vNoise_iniX = 9;vNoise_endX = 51;

    % --- B) - get poinst for vials 128by128
    v1_iniY  = 37;  v1_endY  = 40;      v1_iniX  = 63; v1_endX  = 67;
    
    vNoise_iniY = 108; vNoise_endY = 126;     vNoise_iniX = 102;vNoise_endX = 128;        
    
    % 5.3 - get masks for vials
    Vial_mask(v1_iniX:v1_endX-1,v1_iniY:v1_endY-1,1)      = ones(size(1:v1_endX-v1_iniX,2),size(1:v1_endY-v1_iniY,2));

        % 5.3.2 - SNR 
    imgSE      = permute(imgSE_final,[2 1 4 3]); % struct - [Ny, Nx, #echos, #slices]    
    Vial_noise(vNoise_iniX:vNoise_endX-1,vNoise_iniY:vNoise_endY-1) = ones(size(1:vNoise_endX-vNoise_iniX,2),size(1:vNoise_endY-vNoise_iniY,2));    

            % test SNR
    for echo=1:size(imgSE,3)
        aux_meanNoise(:,:,echo)   = Vial_noise.*imgSE(:,:,echo);
        meanNoise_Intensity(echo) = mean(nonzeros( aux_meanNoise(:,:,echo)));        
        stdVial_T2map_exp(nv)     = std(nonzeros( Vial_mask.*imgSE(:,:,echo)));
    end
    
            % test sigma
    for vi=1:size(Vial_mask,3)
        for echo=1:size(imgSE,3)
            imgSE_norm(:,:,echo)           = imgSE(:,:,echo)./max(max( imgSE(:,:,echo) ));        
            aux_meanVial_norm(:,:,echo,vi) = Vial_mask(:,:,vi).*imgSE_norm(:,:,echo);
            meanVial_norm(echo,vi)         = mean(nonzeros( aux_meanVial_norm(:,:,echo,vi) ));
            sigmaVial_norm(echo,vi)        = std(nonzeros( aux_meanVial_norm(:,:,echo,vi) ));
            SNRVial_norm(echo,vi)          = 1/sigmaVial_norm(echo,vi);
        end
    end
    
    Table_leg    = {'T2 theor.';'SNR echo1';'SNR echo2';'SNR echo3';'SNR echo4';'SNR echo5';'SNR echo6'};
    vial_1      = [581.33;SNRVial_norm(1,1);SNRVial_norm(2,1);SNRVial_norm(3,1);SNRVial_norm(4,1);SNRVial_norm(5,1);SNRVial_norm(6,1)];
                
    T = table(Table_leg,vial_1)  
                
        % ...
    
        % 5.3.3 - Figures
    figure(); imshow(imgSE(:,:,2),[0 max(max(imgSE(:,:,9)))]); %display image 1st echo random slice with good cartilage display
    
    % 5.4 - get T2 maps
    if saveResults == 'True'
        tic 
        % 5.4.1 - initialize
        T2_exp_map = zeros(1,nvials);
            
        for nv = 1:nvials            
            % 5.4.2 - get average per vial
                aux_meanVial(:,:,echo,nv)   = Vial_mask(:,:,nv).*imgSE(:,:,echo);
            
            for echo=1:size(imgSE,3)
                aux_meanVial(:,:,echo,nv)   = Vial_mask(:,:,nv).*imgSE(:,:,echo);
                meanVial_Intensity(echo,nv) = mean(nonzeros( aux_meanVial(:,:,echo,nv)));
                SNR(echo,nv)                = meanVial_Intensity(echo,nv)/std(nonzeros( aux_meanVial(:,:,echo,nv)));
                
            end
            
% %             figure();plot(meanVial_Intensity)
            % 5.4.3 - get mono-exponential fit            
            [T2_exp_map(1,nv) , fit_out{1,nv}] = monoExpfit(t(echoStart:echoEnd(nv)),double(meanVial_Intensity(echoStart:echoEnd(nv),nv)));
            fprintf(['   ----  Mono-exponential fit for Point ', num2str(nv),' / ',num2str(nvials),'   ------ \n'])
        end

        % 5.4.3 - save results
        cd(data_dir)    
        save(['resultsT2_avgBefore_monoexp_maps.mat'],'T2_exp_map')
        if nTE>2
            save(['resultsT2_avgBefore_monoexp_maps.mat'],'T2_exp_map')
        end
    else
        load('resultsT2_avgBefore_monoexp_maps.mat')
    end

    % .. get maxValue ..
    ROI = [1];
    % .. Deviation from theorical ..
    dev_T2vial = (theoric_T2Vial - T2_exp_map)./theoric_T2Vial*100;
    
    % .. get maxValue ..
    figure()
    title(['Echoes used - First:' num2str(echoStart) ' to End:' num2str(echoEnd)])
    yyaxis left    
    plot(ROI, T2_exp_map)
    ylabel('T_2  (ms)')
    ylim([min(T2_exp_map)-10 max(T2_exp_map)+10])
    hold on
    yyaxis right
    ylabel('Deviation from NMR reference (%)')
    plot(ROI, dev_T2vial)
    xlabel('ROI')   
    ylim([-225 5])
    grid on

    
    % ... get table ...
    % print values
    label       = {'Theoric';'T2map_exp';'Std';'SNR';'log10(SNR)'};
    vial_1      = [theoric_T2Vial(1);T2_exp_map(1);stdVial_T2map_exp(1);SNR(2,1);10*log10(SNR(2,1))];

    T = table(label,vial_1)
                
    fprintf('\n\n 5 - Sucessfully finished - avg before - T2 mapping from mono-exponential fit\n\n')
    
end

