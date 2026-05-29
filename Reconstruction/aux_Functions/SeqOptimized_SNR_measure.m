% Get SNR from two signals of the same optimzied sequence
clc
clear all
close all

%% 0 - Setup

clear all
clc
close all

slice = 1;

%% 1 - load
% select directory for recon.
myCD = 'C:\Users\filia\Documents\PhD\Projetos\qMRI\Sequences\MSE\';
cd(myCD)

% DATA 1
    % get dir
str = {'select file'};
data_dir1=uigetdir;
cd(data_dir1)

    % load data
data1 = load('GRAPPA_recon_matlab.mat'); % aux_k - size(Nx, coils, Ny, nslice, nechoes)

    % Get Mask
cd([data_dir1 '\results']) 
mask  = load('roiMSE_T2Phantom_py.mat');
cd(data_dir1)

% DATA 2
    % get dir
str = {'select file'};
data_dir2=uigetdir;
cd(data_dir2)

    % load data
data2 = load('GRAPPA_recon_matlab.mat'); % aux_k - size(Nx, coils, Ny, nslice, nechoes)

%% 2 - Get SNR
% Get signal from data per mask
imgFull1  = permute(data1.img_PI,[2 3 1 4 5]);       % coils, Ny, Nx
imgFull2  = permute(data2.img_PI,[2 3 1 4 5]);       % coils, Ny, Nx
[ncoils, Ny, Nx, nslices, nechos] = size(imgFull1);

% Get images for 1 slice
for eco=1:nechos
    recon_Img_1(:,:,eco) = squeeze(sum(abs(imgFull1(:,:,:,slice,eco).^2),1).^.5);  % LS method for coil combination % size(aux_img_PI) = 18 258 256 3 12
    recon_Img_2(:,:,eco) = squeeze(sum(abs(imgFull2(:,:,:,slice,eco).^2),1).^.5);  % LS method for coil combination % size(aux_img_PI) = 18 258 256 3 12

    diffImage(:,:,eco) = recon_Img_1(:,:,eco) - recon_Img_2(:,:,eco);
    meanImage(:,:,eco) = ( recon_Img_1(:,:,eco) + recon_Img_2(:,:,eco)) / 2;
end

% Plots
maxMean = max(meanImage(:));
% Figure of differences
figure() 
for eco=1:nechos
    subplot(4,4,eco)
    imagesc(squeeze(diffImage(:,:,eco) / maxMean))
    title(['Difference | image Echo ' num2str(eco)])
    colormap gray
    caxis([ min( diffImage(:)/maxMean )  max( diffImage(:)/maxMean )])
end

% Figure of mean
figure()
for eco=1:nechos
    subplot(4,4,eco)
    imagesc(squeeze(meanImage(:,:,eco) / maxMean))
    title(['Mean image | Echo ' num2str(eco)])
    colormap gray
    caxis([ min( meanImage(:)/maxMean )   max( meanImage(:)/maxMean ) ])
end

% Get SNR
nvials = size(mask.T2_mask,3);

for vial=1:nvials
    for eco=1:nechos
        aux_mean = squeeze(meanImage(:,:,eco) .* mask.T2_mask(:,:,vial)');
        aux_diff = squeeze(diffImage(:,:,eco) .* mask.T2_mask(:,:,vial)');

        meanImageVial(eco,vial) =  mean( nonzeros( aux_mean(:) ) );
        diffImageVial(eco,vial) =  mean( nonzeros( aux_diff(:) ) );
    end
end

% SNR=10log10(Psignal/Pnoise)
SNRallEchos_allVials = 10*log10( meanImageVial ./ abs(diffImageVial) );

figure()
imagesc(SNRallEchos_allVials)
hold on
colormap jet
ylabel('Echos')
xlabel('Vials')
title('SNR (db)')
caxis([0   45 ])
