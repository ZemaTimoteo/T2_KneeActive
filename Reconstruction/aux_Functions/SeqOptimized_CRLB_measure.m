% Get CRLB from two signals of the same optimzied sequence
clc
clear all
close all

%% 0 - Setup

clear all
clc
close all

slice    = 1;

%% 1 - load
% select directory for recon.
myCD = 'C:\Users\filia\Documents\PhD\Projetos\qMRI\Sequences\MSE\';
cd(myCD)

% DATA
    % get dir
str = {'select file'};
data_dir1=uigetdir;
cd(data_dir1)

    % load data
data1 = load('GRAPPA_recon_matlab.mat'); % aux_k - size(Nx, coils, Ny, nslice, nechoes)

    % Get Mask
cd([data_dir1 '\results']) 
mask   = load('roiMSE_T2Phantom_py.mat');
nvials = size(mask.T2_mask,3);
cd(data_dir1)


%% 2 - Get AUC
% Get signal from data per mask
imgFull1  = permute(data1.img_PI,[2 3 1 4 5]);       % coils, Ny, Nx
[ncoils, Ny, Nx, nslices, nechos] = size(imgFull1);
totechos  = [1:nechos];

% Get images for 1 slice
for eco=1:nechos
    recon_Img_1(:,:,eco) = squeeze(sum(abs(imgFull1(:,:,:,slice,eco).^2),1).^.5);  % LS method for coil combination % size(aux_img_PI) = 18 258 256 3 12
end


% Figure
% % figure() 
% % for eco=1:nechos
% %     subplot(4,4,eco)
% %     imagesc(squeeze(recon_Img_1(:,:,eco) / maxMValue))
% %     title(['image Echo ' num2str(eco)])
% %     colormap gray
% % %     caxis([0  maxMValue])
% % end

signalEvol = zeros(Nx,Ny,nechos,nvials);
for vial=1:nvials
    for eco=1:nechos
        signalEvol(:,:,eco,vial) = squeeze(recon_Img_1(:,:,eco) .* mask.T2_mask(:,:,vial)');
    end
end

maxMValue = max(signalEvol(:));


% Figure
% % figure() 
% % for eco=1:nechos
% %     subplot(4,4,eco)
% %     imagesc(squeeze(signalEvol(:,:,eco,9) / maxMValue))
% %     title(['image Echo ' num2str(eco)])
% %     colormap gray
% % %     caxis([0  maxMValue])
% % end

% Figure of Signal Evol per vial
figure()
for vial=1:nvials
    subplot(4,4,vial)
    [idxX,idxY]  = find(signalEvol(:,:,1,vial)~=0);
    numberpoints = length(idxX);
    for aux=1:numberpoints
        aux_decay(aux,:,vial) = squeeze(signalEvol(idxX(aux),idxY(aux),:,vial)/maxMValue)';
    end
    plot(squeeze(aux_decay(:,:,vial)'))
    title(['Signal Decay | Vial ' num2str(vial)])
end

figure()
for vial=1:nvials
    subplot(4,4,vial)
    [a,~] = find(aux_decay(:,:,vial)==0); % get zeros elements due to mask
    if a~=0
        mean_decay = mean(aux_decay(1:a(1),:,vial),1);
        std_decay  = std(aux_decay(1:a(1),:,vial),1);
        AUC_vials(:,vial) = trapz(totechos,aux_decay(:,:,vial),2);
    else
        mean_decay = mean(aux_decay(:,:,vial),1);
        std_decay  = std(aux_decay(:,:,vial),1);        
        AUC_vials(:,vial) = trapz(totechos,aux_decay(:,:,vial),2);
    end
    errorbar(totechos,mean_decay,std_decay)
    title(['Signal Decay | Vial ' num2str(vial)])
end

% Figure - AUC
vialsList = {'vial 1','vial 2','vial 3','vial 4','vial 5','vial 6','vial 7','vial 8',...
    'vial 9','vial 10','vial 11','vial 12','vial 13','vial 14',};
AUC_vials(find(AUC_vials==0)) = NaN;
figure()
boxplot(AUC_vials,'Notch','on','Labels',vialsList)
title('Compare Random Data from Different Distributions')
hold on
grid on
ylabel('a.u.')
xlabel('Vials')
title('AUC')
