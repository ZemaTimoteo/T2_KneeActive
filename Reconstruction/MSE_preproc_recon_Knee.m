% MSE_preproc_recon_vF.m function
%     - Complementary function of python code, recon GRAPPA OR LORAKS and creates
%     NIFTI files for denoising and removing Gibbs rings with MRTRIX in
%     Linux.
% 
% by: TTFernandes, July 2024

%% 0 - Setup

clear all
clc
% close all


% toolboxes
cd('D:\Tiago\Trabalho\2021_2025_PhD\Projects\qMRI_Joint\Code\matlabCode\MSE_recon\GRAPPA')
addpath(genpath('D:\Tiago\Trabalho\2021_2025_PhD\Projects\qMRI_Joint\Code\matlabCode\MSE_recon\aux_Functions'));
addpath(genpath('D:\Tiago\Trabalho\2021_2025_PhD\Projects\qMRI_Joint\Code\matlabCode\Toolboxes\UnRing_tool'))
addpath(genpath('D:\Tiago\Trabalho\2021_2025_PhD\Projects\qMRI_Joint\Code\matlabCode\Tutorials\GRAPPA_tutorial'));
addpath(genpath('D:\Tiago\Trabalho\2021_2025_PhD\Projects\qMRI_Joint\Code\matlabCode\Toolboxes\LORAKS_V2.1'))

% test
plotTest     = 'True';
niftTest     = 'Fals';
maskTest     = 'Fals';
ParallelTest = 'LORAK'; % 'GRAPP' or 'SENSE' or 'LORAK'
% filterTest   = 'True';

slice = 1;
echo  = 1;

%% 1 - load
% select directory for recon.
myCD = 'C:\Users\filia\Documents\PhD\Projetos\qMRI\Knee\';
cd(myCD)

% get dir
str = {'select file'};
data_dir=uigetdir;
cd(data_dir)

mkdir([data_dir '\results'])

% load file
if ParallelTest == 'GRAPP'
    load('GRAPPA_for_data_test.mat') % aux_k - size(Nx, coils, Ny, nslice, nechoes)
    aux_k = double(aux_k);
    
    % Plot data
    if plotTest == 'True'
        figure();
        coil = 12;                                          % coil for plot
        imagesc(abs(squeeze(aux_k(:,coil,:,1,1))));
        hold on, title(['Img - K-space'])
    end

elseif ParallelTest == 'LORAK'

    % Load Data
    load('LORAKS_for_data_test.mat') % aux_k - size(Nx, coils, Ny, nslice, nechoes)
    aux_k              = double(aux_k);
    sparseMask_Kspace = double(mask_Kspace);
    
    % Get Full Matrices
    [Nx,ncoils,Ny_real,nsli,ETL] = size(aux_k);
    full_Kspace                     = zeros(Nx, Nx, ncoils,ETL,nsli);  % Get Kspace Full
    fullMask_Kspace                 = zeros(Nx, Nx, ncoils,ETL,nsli);  % Get mask Kspace Full
    perm_aux_k                      = permute(aux_k,[3,1,2,4,5]);

    for sli=1:nsli
        for echo = 1:ETL
            Ny_vector                             = find(sparseMask_Kspace(:,echo)==1);
            full_Kspace(Ny_vector,:,:,echo,sli)   = perm_aux_k(:,:,:,sli,echo);
            fullMask_Kspace(Ny_vector,:,:,echo,:) = 1;
        end
    end

    % Figures
    if plotTest
        figure()
        subplot(121)
        montage(mat2gray(abs(squeeze(full_Kspace(:,:,1,:,1)))))
        title('Full Kspace for all Echoes')        
        subplot(122)
        montage(mat2gray(squeeze(fullMask_Kspace(:,:,1,:,1))))
        title('Full Kspace Mask for all Echoes')
    end

end

%% 2 - Recon Acceleration
tic
%% 2.1 - SENSE
if ParallelTest == 'SENSE'
    
    
 %% 2.2 - GRAPPA
 elseif ParallelTest == 'GRAPP'
    % -------------- 2 - Get data for GRAPPA ----------    
    % 2.1 - Parameters
    sig_red  = permute(aux_k,[2 3 1 4 5]);       % coils, Ny, Nx
    af       = double(R);                    % Acceleration data
%     acs      = sig_red(:, iniKfull+1:endKfull, size(sig_red,3)/2+1-fullLin/2:size(sig_red,3)/2+fullLin/2, :, :); % coils, Ny, Nx, nslice, nechoes    
    acs      = sig_red(:, iniKfull+1:endKfull, :, :, :); % coils, Ny, Nx, nslice, nechoes    
    ksize_x  = 3; % kernel size
    ksize_y  = 2; % kernel size
    
    figure()
    imagesc(abs(squeeze(acs(coil,:,:,1,1))));
    clear aux_k_csm coil_map alias_img           
    
    % ------------- 3 - Apply GRAPPA --------------------------
    % Inputs:
    %   -> sig_red - Reduced data set      (Nc,Ny/r,Nx)  (Nc - number of coils)
    %   -> acs     - Coil sensitivity maps (Nc,Nx,Ny) reduced dimentions
    %   -> af      - Reducion Factor of GRAPPA
    %   -> ksize_x - Kernel size of xx
    %   -> ksize_y - Kernel size of yy
    % Ouputs
    %   -> recon - reconstructed image (Nx,Ny)
    
    % --- 3.1.0 Initialize ---
    [ncoils,orgSizeNx,Nx,nsli,ETL] = size(sig_red);
    Ny                             = Nx;
    sigrecon                       = zeros(ncoils,Nx,Ny,nsli,ETL);
    for ii=1:nsli
        for jj=1:ETL
            % 3.0 - Get kdata
            kdata = sig_red(:,[[1:iniKfull+1] [iniKfull+2:af:endKfull+2] [endKfull+1+2:orgSizeNx]],:,ii,jj); % the plus 1 is because of python index, the +2 is to get new position
                       
            % 3.1 - Apply GRAPPA
            [~,~,sigrecon] = opengrappa_test( ...
                kdata, acs, af, ksize_x, ksize_y ...
                );


            % 3.2 - get recon signal & image
            sigrecon_full = sigrecon;
            
            % 3.3 - Recon data
            out      = ifftshift(ifftshift(ifft(ifft(ifftshift(ifftshift(sigrecon,2),3),[],2),[],3),2),3);
            
            % 3.4 - Export
            aux_kdata_PI(:,:,:,ii,jj) = sigrecon_full;
            aux_img_PI(:,:,:,ii,jj)   = out;

            % 3.5 - Get data from all coils
            recon_Img_full(:,:,ii,jj) = squeeze(sum(abs(aux_img_PI(:,:,:,ii,jj).^2),1).^.5);  % LS method for coil combination % size(aux_img_PI) = 18 258 256 3 12
                        
            % save nift
            if niftTest == 'True'
                v(:,:,ii,jj) = recon_Img_full;

                dir_prefilter = [data_dir '\imagResults_orig'];
                mkdir(dir_prefilter)
                cd(dir_prefilter)
                niftiwrite(v,'imageRecon_prefilter.nii');
                cd(data_dir)
            end
            
            fprintf(['  --  GRAPPA recon for echo ',num2str(jj),' / ',num2str(ETL),'\n\n'])            
        end
    end
    
elseif ParallelTest == 'LORAK'
    % -------------- 2 - Get data for LORAK ----------    
    % 2.1 - Parameters
    accFactor  = double(R);
    rankLORAKS = 45;
    Ny         = Nx;
    sigrecon   = zeros(ncoils,Nx,Ny,nsli,ETL);
    
    nsli = 1; % Avoid recon in all slices
    % ------------- 3 - Apply LORAKS --------------------------
    % Inputs:
    %   -> kdata            - Reduced data set      (Nx, Ny, ncoils,ETL,nsli)  
    %   -> mask             - Mask of full_kspace   (Nx, Ny, ncoils,ETL,nsli)
    %   -> rankLORAKS       - Reducion Factor of LORAKS
    % Ouputs
    %   -> Full_Kspace - reconstructed image (Nx,Ny)
    for sli=1:nsli
        for coil=1:ncoils
            % 3.0 - Get kdata
            kdata = full_Kspace(:, :, coil,:,sli);
            mask  = fullMask_Kspace(:, :, coil,:,sli);

            % 3.1 - Apply LORAKS
            sigrecon(coil,:,:,sli,:) = P_LORAKS( kdata, mask, rankLORAKS );  % Under sampled Data, Mask, Rank
%             sigrecon(coil,:,:,sli,:) = AC_LORAKS( kdata, mask, rankLORAKS );  % Under sampled Data, Mask, Rank

            fprintf(['  --  LORAKS recon slice: ',num2str(sli) ,' / ' ,num2str(nsli) ,' for coil: ',num2str(coil),' / ',num2str(ncoils),'\n\n'])            
            
        end

        for echo=1:ETL
            % 3.2 - get recon signal & image
            sigrecon_full = sigrecon(:,:,:,sli,echo);
            
            % 3.3 - Recon data
            out      = ifftshift(ifftshift(ifft(ifft(ifftshift(ifftshift(sigrecon(:,:,:,sli,echo),2),3),[],2),[],3),2),3);
            
            % 3.4 - Export
            aux_kdata_PI(:,:,:,sli,echo) = sigrecon_full;
            aux_img_PI(:,:,:,sli,echo)   = out;

            % 3.5 - Get data from all coils
            recon_Img_full(:,:,sli,echo) = squeeze(sum(abs(aux_img_PI(:,:,:,sli,echo).^2),1).^.5);  % LS method for coil combination % size(aux_img_PI) = 18 258 256 3 12
% %             recon_Img_full(:,:,ii,jj) = squeeze(abs(sum(aux_img_PI(:,:,:,slice,echo).*conj(aux_img_PI(:,:,:,slice,echo)),ndims(aux_img_PI(:,:,:,slice,echo)))).^(1/2));
                        
            % save nift
            if niftTest == 'True'
                v(:,:,sli,echo) = recon_Img_full;

                dir_prefilter = [data_dir '\imagResults_orig'];
                mkdir(dir_prefilter)
                cd(dir_prefilter)
                niftiwrite(v,'imageRecon_prefilter.nii');
                cd(data_dir)
            end
            
            fprintf(['  --  LORAKS recon for echo ',num2str(echo),' / ',num2str(ETL),'\n\n'])            
        end
    end
    
end

toc

%% ----------- 4 - Get Figures -----------------------------
% close all
if plotTest == 'True'
    cd([data_dir '\results'])
    % 4.1 - Kspace
    figure() 
    kspace_plot = squeeze(abs(aux_kdata_PI(1,:,:,slice,echo)));
    imagesc(kspace_plot)
    hold on, title(['K space'])
    colormap gray

    % 4.2 - Recon Images
    fig1 = figure();
    for sli =1:nsli
        plotImageRecon = squeeze(sum(abs(aux_img_PI(:,:,:,sli,echo).^2),1).^.5);  % LS method for coil combination
        subplot(1,nsli,sli)
        imagesc(plotImageRecon,[0 .4])
        hold on, title(['Slice = ',num2str(sli)])
        colormap gray
        caxis([ 0   max( max(plotImageRecon) ) ])
    end    
    saveas(fig1,'Recon_image_&_echos.png')



    % 4.3 - Montage per echo
    fig2 = figure();
    subplot(121)
    h = montage(mat2gray(squeeze(abs(recon_Img_full(:,:,1,:)))));
    title(['All Echos - Recon ' ParallelTest])
    subplot(122)
    imagesc(recon_Img_full(:,:,1,1))
    colormap gray
    title(['1st Echos - 1st Slice Recon ' ParallelTest])
    saveas(fig2,'Recon_1slice.png')

end
    

%% -------------- 4.5 - Get mask  ------------------------------    
if maskTest == 'True'
    
    nvials = 14;
    for jj=1:nvials
        % --  Circular Vial mask     
        T2_mask_m = zeros(size(recon_Img_full(:,:,1)));
        figure('Name','Results','WindowState','maximized'); imshow(recon_Img_full(:,:,1,3),[]);           % display image 1st echo random slice with good cartilage display
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

        % --  Noise mask
% %             figure; imshow(recon_Img_full(:,:,1),[]); % display image 1st echo random slice with good cartilage display
% %             title('Create Mask for Preprocessing')
% %             Noise_mask_m  = roipoly;                      % select roi with cartilage region only
% %             Noise_mask    = double(Noise_mask_m);       
    end
    
    % --  save mask
    cd([data_dir '\results'])
    name_mask = 'roiMSE_T2Phantom_py.mat';
    save(name_mask,'T2_mask')
    
    
end
            

%% ------------ 5 - Save kdata and img ----------------------
cd(data_dir)
kdata_PI = permute(aux_kdata_PI,[2 1 3 4 5]); % Nx, Ny, coils, nslice, #echoes
img_PI   = permute(aux_img_PI,[2 1 3 4 5]);   % Nx, Ny, coils, nslice, #echoes

if ParallelTest == 'SENSE'

elseif ParallelTest == 'GRAPP'
    save('GRAPPA_recon_matlab.mat','kdata_PI','img_PI')
    filename = 'reconGRAPPA_imag';
    fprintf(['  --  GRAPPA recon Performed and Saved  -- \n\n'])
elseif ParallelTest == 'LORAK'
    save('LORAKS_recon_matlab.mat','kdata_PI','img_PI')
    filename = 'reconLORAKS_imag';
    fprintf(['  --  LORAKS recon Performed and Saved  -- \n\n'])
end
    
    %% -------------- 6 - Filter image ------------------------------
% % %     if filterTest == 'True'
% % %         % 5.0 - normal
% % %         figure()
% % %         subplot(221)
% % %         imagesc(recon_Img_full,[0 .4])
% % %         colormap gray
% % %         hold on, title(['Original'])
% % %         caxis([ 0   max( max(recon_Img_full) ) ])
% % %         
% % %         % 5.1 - Gaussian filter  Image Space
% % %         sigma = 2;
% % %         testIMG_1 = imgaussfilt(recon_Img_full,sigma);
% % %         subplot(222)
% % %         imagesc(testIMG_1,[0 .4])
% % %         colormap gray
% % %         hold on, title(['Image space Gaussian sigma=',num2str(sigma)])
% % %         caxis([ 0   max( max(recon_Img_full) ) ])
% % %         
% % %         % 5.2 - Gaussian filter  K-space
% % %         sigma = 0.8;
% % %         testIMG_1 = imgaussfilt(recon_Img_full,sigma);
% % %         subplot(222)
% % %         imagesc(testIMG_1,[0 .4])
% % %         colormap gray
% % %         hold on, title(['Image space Gaussian sigma=',num2str(sigma)])
% % %         caxis([ 0   max( max(recon_Img_full) ) ])
% % %         
% % %         
% % %         % 5.2 --- Ringing artefact removal ---
% % %         % 5.2.1 - HAMMING FILTER
% % %         for i=1:ncoils
% % %             kspace_plot = squeeze(abs(aux_kdata_PI(i,:,:,1,1)));
% % %             
% % %             Ndim = ndims(kspace_plot);
% % %             prec = class( kspace_plot );
% % %             filter = cast( 1, prec );
% % %             sz = size( kspace_plot );
% % %             for d = 1:Ndim % Excluding channel dimension
% % %                 reshRule = ones(1,Ndim);    % how the filter will be reshaped
% % %                 filt_1Dime = hamming( sz(d), 'periodic');
% % %                 reshRule(d) = sz(d);
% % %                 filt_1Dime = reshape( filt_1Dime, reshRule );
% % %                 filter = bsxfun( @times, filter, filt_1Dime ); % Build up mutli dimensional fiter
% % %             end
% % %             kspFiltered(i,:,:) = bsxfun( @times, filter, kspace_plot); % Apply fiter on kspace data
% % %             clear kspace_plot
% % %         end
% % %         %Reconstruction and channel combine - After ringing artefact removal
% % %         aux_recon_Img_filtHAMM  = ifftshift(ifftshift(ifft(ifft(ifftshift(ifftshift(kspFiltered,2),3),[],2),[],3),2),3);
% % %         recon_Img_full_filtHAMM = squeeze(sum(abs(aux_recon_Img_filtHAMM.^2),1).^.5);  % LS method for coil combination
% % %         subplot(223)
% % %         imagesc(recon_Img_full_filtHAMM,[0 .4])
% % %         hold on, title(['Hamming Filter'])
% % %         colormap gray
% % %         caxis([ 0   max( max(recon_Img_full) ) ])
% % %         
% % %         % 5.2.2 - Kaiser FILTER
% % %         % Filter Definition
% % %         beta = 1; % roughly equivalent to Hann
% % %         w    = kaiser(size(aux_kdata_PI(1,:,:,1,1),2),beta);
% % %         
% % %         % sHR: HR reference image
% % %         n     = size(out_full,2);
% % %         for i=1:ncoils
% % %             kspace_plot(i,:,:)        = fftshift(fft(squeeze(out_full(i,:,:)),n,1),1)/n; % HR image in k-space
% % %             kspFiltered_kaiser(i,:,:) = w .* squeeze(kspace_plot(i,:,:));
% % %             aux_img(i,:,:)               = real(ifft(ifftshift(squeeze(kspFiltered_kaiser(i,:,:)),1),n,1))*n;
% % %         end
% % %         recon_Img_full_filtHKAI = squeeze(sum(abs(aux_img.^2),1).^.5);  % LS method for coil combination
% % %         
% % %         subplot(224)
% % %         imagesc(recon_Img_full_filtHKAI,[0 .4])
% % %         hold on, title(['Kaiser Filter'])
% % %         colormap gray
% % %         caxis([ 0   max( max(recon_Img_full) ) ])
% % %     end
% % %     
% % %     % 5.4 - Filter with unring (tool for removal of the Gibbs ringing artefact - Kellner, et al. 2016 MRM)
% % %     % 5.4.1 - Parameters
% % %     % % % kdata_Test   = aux_kdata_PI;
% % %     % % % v      = squeeze(kdata_Test(1,:,:,1,1));                  % input volume
% % %     % % % nsh    = 20;                % discretization of subpixel spaceing (default 20)
% % %     % % % minW   = 1;                 % left border of window used for TV computation (default 1)
% % %     % % % maxW   = 3;                 % right border of window used for TV computation (default 3)
% % %     % % % params = [nsh minW maxW];
% % %     % % %
% % %     % % % v = unring(v,params);
% % %     
    
    %% ------------ 7 - Check SNR values ----------------------    

% SNR test 304 d12_07_22    
% % % % cd('C:\Users\filia\Documents\PhD\Projetos\qMRI\Sequences\MSE\test304_MSE_d12_07_22v2\results')
% % % % load('reconData.mat')
% % % % load('roiMSE_T2Phantom_py.mat')    
% % % % clear data_Vial MeanVial stanDeviation MeanVial_norm stanDeviation_norm SNR_vial_norm SNR_vial norm_all_recon_Img_full
% % % % all_recon_Img_full = squeeze(image);
% % % % data_Vial     = zeros(Nx,Ny/2,ETL,size(T2_mask,3));
% % % % SNR_vial      = zeros(ETL,size(T2_mask,3));
% % % % SNR_vial_norm = zeros(ETL,size(T2_mask,3));
% % % % for ii=1:size(T2_mask,3)
% % % %     for jj = 1:size(all_recon_Img_full,3)
% % % %         data_Vial(:,:,jj,ii) = all_recon_Img_full(:,:,jj) .* T2_mask(:,:,ii) ;        
% % % %         stanDeviation(jj,ii) = std(nonzeros(data_Vial(:,:,jj,ii)));
% % % %         MeanVial(jj,ii)      = mean(nonzeros(data_Vial(:,:,jj,ii)));
% % % %         SNR_vial(jj,ii)      = MeanVial(jj,ii)/stanDeviation(jj,ii);
% % % %         
% % % %         norm_all_recon_Img_full(:,:,jj) = all_recon_Img_full(:,:,jj)./max(max(all_recon_Img_full(:,:,jj)));        
% % % %         data_Vial_norm(:,:,jj,ii) = norm_all_recon_Img_full(:,:,jj) .* T2_mask(:,:,ii) ;                        
% % % %         stanDeviation_norm(jj,ii) = std(nonzeros(data_Vial_norm(:,:,jj,ii)));
% % % %         MeanVial_norm(jj,ii)      = mean(nonzeros(data_Vial_norm(:,:,jj,ii)));
% % % %         SNR_vial_norm(jj,ii)      = MeanVial_norm(jj,ii)/stanDeviation_norm(jj,ii);
% % % %     end
% % % % end
% % % % 
% % % %     % print values
% % % %     T2map       = {'T2 theor.';'SNR echo1';'SNR echo2';'SNR echo3';'SNR echo4';'SNR echo5';'SNR echo6'};
% % % %     vial_1      = [581.33;SNR_vial_norm(1,1);SNR_vial_norm(2,1);SNR_vial_norm(3,1);SNR_vial_norm(4,1);SNR_vial_norm(5,1);SNR_vial_norm(6,1)];
% % % %     vial_2      = [403.5 ;SNR_vial_norm(1,2);SNR_vial_norm(2,2);SNR_vial_norm(3,2);SNR_vial_norm(4,2);SNR_vial_norm(5,2);SNR_vial_norm(6,2)];
% % % %     vial_3      = [278.1 ;SNR_vial_norm(1,3);SNR_vial_norm(2,3);SNR_vial_norm(3,3);SNR_vial_norm(4,3);SNR_vial_norm(5,3);SNR_vial_norm(6,3)];
% % % %     vial_4      = [190.94;SNR_vial_norm(1,4);SNR_vial_norm(2,4);SNR_vial_norm(3,4);SNR_vial_norm(4,4);SNR_vial_norm(5,4);SNR_vial_norm(6,4)];
% % % %     vial_5      = [133.27;SNR_vial_norm(1,5);SNR_vial_norm(2,5);SNR_vial_norm(3,5);SNR_vial_norm(4,5);SNR_vial_norm(5,5);SNR_vial_norm(6,5)];
% % % %     vial_6      = [96.89 ;SNR_vial_norm(1,6);SNR_vial_norm(2,6);SNR_vial_norm(3,6);SNR_vial_norm(4,6);SNR_vial_norm(5,6);SNR_vial_norm(6,6)];
% % % %     vial_7      = [64.07 ;SNR_vial_norm(1,7);SNR_vial_norm(2,7);SNR_vial_norm(3,7);SNR_vial_norm(4,7);SNR_vial_norm(5,7);SNR_vial_norm(6,7)];
% % % %     vial_8      = [46.42 ;SNR_vial_norm(1,8);SNR_vial_norm(2,8);SNR_vial_norm(3,8);SNR_vial_norm(4,8);SNR_vial_norm(5,8);SNR_vial_norm(6,8)];
% % % %     vial_9      = [31.97 ;SNR_vial_norm(1,9);SNR_vial_norm(2,9);SNR_vial_norm(3,9);SNR_vial_norm(4,9);SNR_vial_norm(5,9);SNR_vial_norm(6,9)];
% % % %     vial_10     = [22.56 ;SNR_vial_norm(1,10);SNR_vial_norm(2,10);SNR_vial_norm(3,10);SNR_vial_norm(4,10);SNR_vial_norm(5,10);SNR_vial_norm(6,10)];
% % % %     vial_11     = [15.813;SNR_vial_norm(1,11);SNR_vial_norm(2,11);SNR_vial_norm(3,11);SNR_vial_norm(4,11);SNR_vial_norm(5,11);SNR_vial_norm(6,11)];
% % % %     vial_12     = [11.237;SNR_vial_norm(1,12);SNR_vial_norm(2,12);SNR_vial_norm(3,12);SNR_vial_norm(4,12);SNR_vial_norm(5,12);SNR_vial_norm(6,12)];
% % % %     vial_13     = [7.911 ;SNR_vial_norm(1,13);SNR_vial_norm(2,13);SNR_vial_norm(3,13);SNR_vial_norm(4,13);SNR_vial_norm(5,13);SNR_vial_norm(6,13)];
% % % %     vial_14     = [5.592 ;SNR_vial_norm(1,14);SNR_vial_norm(2,14);SNR_vial_norm(3,14);SNR_vial_norm(4,14);SNR_vial_norm(5,14);SNR_vial_norm(6,14)];
% % % %                 
% % % %     T = table(T2map,vial_1,vial_2,vial_3,vial_4,vial_5,vial_6,vial_7,vial_8,...
% % % %                     vial_9,vial_10,vial_11,vial_12,vial_13,vial_14)  
% % % %                 
% % % %                 
% % % % % ... 
% % % % 
% % % % clear all_recon_Img_full norm_all_recon_Img_full
% % % % for ii = 1:ETL    
% % % %     all_recon_Img_full(:,:,ii)      = squeeze(sum(abs(aux_img_PI(:,:,:,1,ii).^2),1).^.5);  % LS method for coil combination
% % % %     norm_all_recon_Img_full(:,:,ii) = all_recon_Img_full(:,:,ii)./max(max(all_recon_Img_full(:,:,ii)));
% % % % end
% % % % 
% % % % % SNR = médie(vial) / std(vial)
% % % % %  6.1 -> All vial where vial is hommogeneous
% % % % load([data_dir '/results/roiMSE_T2Phantom_py.mat'])
% % % % data_Vial     = zeros(Nx,Ny,ETL,size(T2_mask,3));
% % % % SNR_vial_norm = zeros(ETL,size(T2_mask,3));
% % % % for ii=1:size(T2_mask,3)
% % % %     for jj = 1:size(all_recon_Img_full,3)  
% % % %         data_Vial(:,:,jj,ii) = all_recon_Img_full(:,:,jj) .* T2_mask(:,:,ii) ;        
% % % %         stanDeviation(jj,ii) = std(nonzeros(data_Vial(:,:,jj,ii)));
% % % %         MeanVial(jj,ii)      = mean(nonzeros(data_Vial(:,:,jj,ii)));
% % % %         SNR_vial_norm(jj,ii) = MeanVial(jj,ii)/stanDeviation(jj,ii);
% % % %     end
% % % % end
% % % % 
% % % %     % print values
% % % %     T2map       = {'T2 theor.';'SNR echo1';'SNR echo2';'SNR echo3';'SNR echo4';'SNR echo5';'SNR echo6'};
% % % %     vial_1      = [581.33;SNR_vial_norm(1,1);SNR_vial_norm(2,1);SNR_vial_norm(3,1);SNR_vial_norm(4,1);SNR_vial_norm(5,1);SNR_vial_norm(6,1)];
% % % %     vial_2      = [403.5 ;SNR_vial_norm(1,2);SNR_vial_norm(2,2);SNR_vial_norm(3,2);SNR_vial_norm(4,2);SNR_vial_norm(5,2);SNR_vial_norm(6,2)];
% % % %     vial_3      = [278.1 ;SNR_vial_norm(1,3);SNR_vial_norm(2,3);SNR_vial_norm(3,3);SNR_vial_norm(4,3);SNR_vial_norm(5,3);SNR_vial_norm(6,3)];
% % % %     vial_4      = [190.94;SNR_vial_norm(1,4);SNR_vial_norm(2,4);SNR_vial_norm(3,4);SNR_vial_norm(4,4);SNR_vial_norm(5,4);SNR_vial_norm(6,4)];
% % % %     vial_5      = [133.27;SNR_vial_norm(1,5);SNR_vial_norm(2,5);SNR_vial_norm(3,5);SNR_vial_norm(4,5);SNR_vial_norm(5,5);SNR_vial_norm(6,5)];
% % % %     vial_6      = [96.89 ;SNR_vial_norm(1,6);SNR_vial_norm(2,6);SNR_vial_norm(3,6);SNR_vial_norm(4,6);SNR_vial_norm(5,6);SNR_vial_norm(6,6)];
% % % %     vial_7      = [64.07 ;SNR_vial_norm(1,7);SNR_vial_norm(2,7);SNR_vial_norm(3,7);SNR_vial_norm(4,7);SNR_vial_norm(5,7);SNR_vial_norm(6,7)];
% % % %     vial_8      = [46.42 ;SNR_vial_norm(1,8);SNR_vial_norm(2,8);SNR_vial_norm(3,8);SNR_vial_norm(4,8);SNR_vial_norm(5,8);SNR_vial_norm(6,8)];
% % % %     vial_9      = [31.97 ;SNR_vial_norm(1,9);SNR_vial_norm(2,9);SNR_vial_norm(3,9);SNR_vial_norm(4,9);SNR_vial_norm(5,9);SNR_vial_norm(6,9)];
% % % %     vial_10     = [22.56 ;SNR_vial_norm(1,10);SNR_vial_norm(2,10);SNR_vial_norm(3,10);SNR_vial_norm(4,10);SNR_vial_norm(5,10);SNR_vial_norm(6,10)];
% % % %     vial_11     = [15.813;SNR_vial_norm(1,11);SNR_vial_norm(2,11);SNR_vial_norm(3,11);SNR_vial_norm(4,11);SNR_vial_norm(5,11);SNR_vial_norm(6,11)];
% % % %     vial_12     = [11.237;SNR_vial_norm(1,12);SNR_vial_norm(2,12);SNR_vial_norm(3,12);SNR_vial_norm(4,12);SNR_vial_norm(5,12);SNR_vial_norm(6,12)];
% % % %     vial_13     = [7.911 ;SNR_vial_norm(1,13);SNR_vial_norm(2,13);SNR_vial_norm(3,13);SNR_vial_norm(4,13);SNR_vial_norm(5,13);SNR_vial_norm(6,13)];
% % % %     vial_14     = [5.592 ;SNR_vial_norm(1,14);SNR_vial_norm(2,14);SNR_vial_norm(3,14);SNR_vial_norm(4,14);SNR_vial_norm(5,14);SNR_vial_norm(6,14)];
% % % %                 
% % % %     T = table(T2map,vial_1,vial_2,vial_3,vial_4,vial_5,vial_6,vial_7,vial_8,...
% % % %                     vial_9,vial_10,vial_11,vial_12,vial_13,vial_14)    
% % % % 
% % % % %  6.2 -> SNR pixel to pixel
% % % % SNR_pixel = mean(all_recon_Img_full,3) / squeeze(std(permute(all_recon_Img_full,[3 1 2])));
% % % % SNR_pixel_norm = mean(norm_all_recon_Img_full,3) / squeeze(std(permute(norm_all_recon_Img_full,[3 1 2])));
% % % % 
% % % % %  6.3 -> Figures
% % % % % figure()
% % % % % % imagesc(squeeze(all_recon_Img_full(:,:,2)),[])
% % % % % imagesc(SNR_pixel_norm)
% % % % % imagesc(mean(all_recon_Img_full,3),[0 3])
% % % % % hold on, title(['SNR per pixel over echoes'])
% % % % % colormap gray
% % % % % caxis([ 0   max( max(SNR_pixel) ) ])
% % % % 

%% ------------ 8 - Save Nifti images and mask ---------------

if niftTest == 'True'
    mkdir([data_dir '\imagResults_preproc'])
    % Combine coils
    for eco=1:size(aux_img_PI,5)
        recon_Img_full(:,:,eco) = squeeze(sum(abs(aux_img_PI(:,:,:,slice,echo).^2),1).^.5);  % LS method for coil combination
    end
    
    % Save images as Nifti
    cd([data_dir '\imagResults_orig'] )
    niftiwrite(recon_Img_full,filename)
    
    % Get mask & save as Nifti
    if maskTest == 'True'
        % save nifti
        filename = 'mask_imag';
        niftiwrite(T2_mask,filename)
        fprintf(['  --   Mask as Nifti Performed and Saved  -- \n\n'])            
    end
    cd(data_dir)
end

%% 9 - END
load('sequence_info.mat')
fprintf(['  --  MSE_preproc_recon_vf.m - for test' num2str(test) ' Finnished!  -- \n\n'])

load handel
sound(y,Fs)
