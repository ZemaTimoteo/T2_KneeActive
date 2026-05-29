% Read T2 maps for Repeatability study
clear all 
clc

%% 1 - LOAD Data% select directory for recon.
myCD = 'C:\Users\filia\Documents\PhD\Projetos\qMRI\Sequences\UTE\';
cd(myCD)

% get dir
str = {'select file'};
data_dir=uigetdir;
cd(data_dir)
data_listing = dir(data_dir);

% load file
im_data_file = [data_dir filesep data_listing(5).name];   % Dicom Data
load(im_data_file);

%% 2 - Plot
figure()
for ii=1:size(imspace,3) %where N is the number of images
    subplot(5,4,ii)
    imshow(squeeze(imspace(:,:,ii)),[])
end
% hold on; 
% colormap jet;
% caxis([0 1000])
% colorbar

%% 3 - Create Video
video = VideoWriter('Echoes.avi'); %create the video object
open(video); %open the file for writing
for ii=1:size(imspace,3) %where N is the number of images
    % convert the image to a frame
    frame = im2frame(abs(imspace(:,:,ii)));
    
    writeVideo(video,frame); %write the image to file
end
close(video); %close the file