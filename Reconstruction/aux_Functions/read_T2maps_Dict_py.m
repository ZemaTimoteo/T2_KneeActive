


% Read T2 maps for Repeatability study
clear all 
clc

%% 1 - LOAD Data% select directory for recon.
myCD = 'C:\Users\filia\Documents\PhD\Projetos\qMRI\Sequences\MSE\';
cd(myCD)

% get dir
str = {'select file'};
data_dir=uigetdir;
cd(data_dir)
data_listing = dir(data_dir);

% load file
im_data_file = [data_dir filesep data_listing(6).name];   % Dicom Data
T2map        = load(im_data_file);

rgb = T2map.T2_dict_map(:,:,1,1);
%% 2 - Plot
figure()
imshow(squeeze(T2map.T2_dict_map(:,:,1,1)),[])
hold on; colormap jet;
caxis([0 1000])
colorbar


%% 3 - 
% 3.1 - get diamiter
d        = drawline;
pos      = d.Position;
diffPos  = diff(pos);
diameter = hypot(diffPos(1),diffPos(2))

% 3.2 - Initial Attempt to Find Circles
gray_image = im2gray(rgb);
figure();imshow(gray_image)
[centers,radii] = imfindcircles(rgb,[20 25],'ObjectPolarity','dark')

% 3.3 - Increase Detection Sensitivity
[centers,radii] = imfindcircles(rgb,[20 25],'ObjectPolarity','dark', ...
    'Sensitivity',0.988)

% 3.4 - draw circles
figure(); imshow(rgb)
h = viscircles(centers,radii);
