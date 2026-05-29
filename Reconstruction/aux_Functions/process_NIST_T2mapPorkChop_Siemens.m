% Read T2 maps for Repeatability study from Siemens Dicom files and export
% maps
%
% @TTFernandes, November, 2022

clear all 
clc

plotTest = 'Fals';
maskTest = 'True';
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
vT2dic_name = ['vT2_map_ data_' data_dir(end-10:end-7)];

% Load 1st Echo Imaging
dir_1ech          = data_dir;
data_1ech_listing = dir(dir_1ech);
im_1ech_file      = [dir_1ech filesep data_1ech_listing(3).name];   % Dicom Data
recon_Img_1ech    = double((dicomread(im_1ech_file)));
infoDicom_1ech    = (dicominfo(im_1ech_file));


clear data_listing
%% 2 - Plot
if plotTest == 'True'
    figure()
    imshow(recon_Img_full)
    hold on; colormap jet;
    caxis([0 400])
    colorbar
end

%% 3 - Get Mask
nvials =1 ;
if maskTest == 'True'
    
    % --  Get mask
    figure; imshow(recon_Img_full(:,:,1),[]); % display image 1st echo random slice with good cartilage display
    title('Create Mask for Preprocessing')
    T2_mask_m = roipoly;                  % select roi with cartilage region only
    T2_mask   = double(T2_mask_m);   
    
    % --  save mask
    cd([data_dir '\results'])
    name_mask = 'roiMSE_T2Phantom_py.mat';
    save(name_mask,'T2_mask')
    
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
cmap = [0  80];

% T2 dict
FigPorkchop = figure()
ax1=axes;
imshow(recon_Img_full(:,:,1),[]);axis image;    % display image 1st echo random slice with good cartilage display
% title(['vender T2 map ' T2dic_name])
colormap(ax1,'gray')
ax1.Visible = 'off';
hold on

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

colorbar('southoutside','fontsize',20)

cd(dir_results)
saveas(FigPorkchop,['PorkChopSiemensfig.png'])
cd(data_dir)


%% 6 - Save Results
cd(dir_results)
save(vT2dic_name,'vT2_dict_map')

% -- Save in common folder
cd('C:\Users\filia\Documents\PhD\Projetos\qMRI\Sequences\MSE\RepeatabilityStudy_ISMRM23\Results_maps')
save(vT2dic_name,'vT2_dict_map')

