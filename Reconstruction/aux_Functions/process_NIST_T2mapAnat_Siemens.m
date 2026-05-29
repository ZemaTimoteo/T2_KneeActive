% Read T2 maps for Repeatability study from Siemens Dicom files and export
% maps
%
% @TTFernandes, November, 2022

clear all 
clc

plotTest = 'Fals';
maskTest = 'True';
addpath(genpath('D:\Tiago\Trabalho\2021_2025_PhD\Projects\qMRI_Joint\Code\matlabCode\MSE_recon\aux_Functions'))
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
im_data_file   = [data_dir filesep data_listing(3).name];   % Dicom Data
recon_Img_full = (dicomread(im_data_file));
infoDicom      = (dicominfo(im_data_file));

% Get name
vT2dic_name = ['vT2_map_ data_' data_dir(end-18:end-15)];

% Load 1st Echo Imaging
% % dir_1ech          = [data_dir(1:end-6) 'dicoms'];
% % data_1ech_listing = dir(dir_1ech);
% % im_1ech_file      = [dir_1ech filesep data_1ech_listing(3).name];   % Dicom Data
% % recon_Img_1ech    = double((dicomread(im_1ech_file)));
% % infoDicom_1ech    = (dicominfo(im_1ech_file));


clear data_listing
%% 2 - Plot
if plotTest == 'True'
    T2mSiemens = figure('WindowState','maximized');          % display image 1st echo random slice with good cartilage display
    imshow(recon_Img_full)
    hold on; colormap jet;
    caxis([0 850])
    colorbar('fontsize',20)
    cd(dir_results)
    saveas(T2mSiemens,['T2mSiemens.png'])
    cd(data_dir)
end

%% 3 - Get Mask

if maskTest == 'True'
    nvials = 14;
    for jj=1:nvials
        % --  Circular Vial mask
        T2_mask_m = zeros(size(recon_Img_full));
        figure('Name','Results','WindowState','maximized'); imshow(recon_Img_full,[]);           % display image 1st echo random slice with good cartilage display
        title(['Create Mask for Preprocessing vial=', num2str(jj)])
        colormap jet;caxis([0 400]);
        
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


%% 4 - Get T2 values p/ mask
vT2_dict_map        = zeros(1,size(T2_mask,1)*size(T2_mask,2));
resh_T2mask         = reshape(T2_mask,1,size(T2_mask,1)*size(T2_mask,1),nvials);
resh_recon_Img_full = reshape(recon_Img_full,1,size(recon_Img_full,1)*size(recon_Img_full,1));

for jj=1:nvials
    aux_vT2_dict_map(:,:,jj) = resh_T2mask(1,:,jj).*double(resh_recon_Img_full);
end
vT2_dict_map = reshape(aux_vT2_dict_map,size(T2_mask,1),size(T2_mask,2),nvials);

%% 5 - Plot Masks
cmap = [0  650];

% T2 dict
figure()
ax1=axes;
imshow(recon_Img_full, [], 'Colormap', gray(256));
axis image;    % display image 1st echo random slice with good cartilage display
% title(['vender T2 map ' T2dic_name])
ax1.Visible = 'off';
hold on
caxis(cmap)

for jj=1:nvials
    alphadata = T2_mask(:,:,jj);
    
    ax2=axes;
    imagesc(vT2_dict_map(:,:,jj),'AlphaData',alphadata); axis image;
    colormap(ax2,'jet');
    caxis(ax2,cmap);
    ax2.Visible = 'off';
    linkprop([ax1 ax2],'Position');
    hold on
end








%% 6 - Save Results
cd(dir_results)
save(vT2dic_name,'vT2_dict_map')

% -- Save in common folder
cd('C:\Users\filia\Documents\PhD\Projetos\qMRI\Sequences\MSE\RepeatabilityStudy_ISMRM23\Results_maps')
save(vT2dic_name,'vT2_dict_map')

fprintf('\n\n -------- Run Sucessfully finished - All code ---------\n\n')
