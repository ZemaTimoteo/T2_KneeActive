clear all
clc
close all

%% 1 - Create data
% 1.1 - Load data
load cmap_head;

cmap = cmap_head;    % coil sensitivity map
img  = phantom(128); % image
nc   = size(cmap,1); % number of coils

% 1.2 - Get Coil Image
for k=1:nc
    coil_img(k,:,:) = squeeze(cmap(k,:,:)).*img; % cmap by image 
end

% (1.3) EXTRA - Addnoise to coil_img
img_coil_noise=addnoise(coil_img,0.001);

% 1.4 - Get Images
figure();
subplot(121);imagesc(abs(img));
subplot(122);imagesc(abs(squeeze(img_coil_noise(1,:,:))));


%% 2 - Process data
% 2.1 - Get k-space
sig_coil_noise = fftshift(fftshift(fft(fft(fftshift(fftshift(img_coil_noise,2),3),[],2),[],3),2),3);

af = 4;  % aceleration data
sx = 16; % autocalibration lines (ACS)x data size
sy = 16; % autocalibration lines (ACS)y data size

% 2.2 - Get Image subsampled
sig_red = sig_coil_noise(:,1:af:end,:); % Reduced sig_coil_size

% 2.3 - Get ACS autocalibration lines
acs = sig_coil_noise(:,65-sy:64+sy,65-sx:64+sx);  % HOW TO GET ACS? Is this a subsampled of k-space=

% 2.4 - Get Figure
figure();
subplot(131);imagesc(squeeze(abs(sig_coil_noise(1,:,:))))
subplot(132);imagesc(squeeze(abs(sig_red(1,:,:))))
subplot(133);imagesc(squeeze(abs(acs(1,:,:))))

%% 3 - GRAPPA
% 3.0 - Parameters
ksize_x1 = 9; ksize_y1 = 8; % kernel size
ksize_x2 = 7; ksize_y2 = 6; % kernel size
ksize_x3 = 5; ksize_y3 = 4; % kernel size
ksize_x4 = 3; ksize_y4 = 2; % kernel size
af       = 4;               % Acceleration Factor

% 3.1 - Apply GRAPPA
out1 = opengrappa(sig_red,acs,af,ksize_x1,ksize_y1);
out2 = opengrappa(sig_red,acs,af,ksize_x2,ksize_y2);
out3 = opengrappa(sig_red,acs,af,ksize_x3,ksize_y3);
out4 = opengrappa(sig_red,acs,af,ksize_x4,ksize_y4);

% 3.2 - Get Figures
figure()
subplot(2,2,1)
imagesc(squeeze(sum(abs(out1).^2,1).^.5),[0 .4])
hold on, title(['kernel x,y = (',num2str(ksize_x1),',',num2str(ksize_y1),')'])
subplot(2,2,2)
imagesc(squeeze(sum(abs(out2).^2,1).^.5),[0 .4])
hold on, title(['kernel x,y = (',num2str(ksize_x2),',',num2str(ksize_y2),')'])
subplot(2,2,3)
imagesc(squeeze(sum(abs(out3).^2,1).^.5),[0 .4])
hold on, title(['kernel x,y = (',num2str(ksize_x3),',',num2str(ksize_y3),')'])
subplot(2,2,4)
imagesc(squeeze(sum(abs(out4).^2,1).^.5),[0 .4])
hold on, title(['kernel x,y = (',num2str(ksize_x4),',',num2str(ksize_y4),')'])
colormap gray