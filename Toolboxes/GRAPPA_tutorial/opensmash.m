function [img_out, sig_out] = opensmash(sig_in,cmap,af)

%   [img_out, sig_out] = opensmash(sig_in,cmap,af);
%
%   SMASH reconstruction code
%       
%   IN:         sig_in             reduced k-space Signal (#coils, Ky./af, Kx)
%               cmap               coil sensitivity maps    (#coils, Ny,Nx)
%               af                 Acceleration factor      (integer)
%
%   OUT:        img_out            Reconstructed Image     (Ny, Nx)    
%               img_out            Reconstructed Signal    (Ny, Nx)
%   
%   20.05.2008  Felix Breuer
% 


[nc_cmap,ny,nx] = size(cmap);
[nc_red,ky_red,kx] = size(sig_in);

if nc_cmap~=nc_red
    disp('Error! The number of coils has to be the same for both inputs!')
    return;
else
    nc = nc_red;
end

if ky_red*af~=ny && ky_red~=ny || kx~=nx
    disp('Error! Size of reduced k-space signal is not compatible with size of coil maps!')
    return;
end

cmap(isinf(cmap))=0;        % Just check to make sure there are no bad entries in the coil maps
cmap(isnan(cmap))=0;


sig_hb = ifftshift(ifft(ifftshift(sig_in,3),[],3),3);           % fft k-space signal (ky_red,kx) into hybrid signal (ky_red,nx) 


dky = 2*pi/ny;              % distance between k-space lines 

y = [-(ny/2)+1:ny/2];       

for m = 1:af,
    harm(m,:) = exp(-i*dky*(m-1)*y);                   % harm     -> spatial harmonics of order m-1  (af,Ny)      
end


for x=1:nx,
    cmap_tmp = squeeze(cmap(:,:,x));                   % cmap_tmp -> sensitivity profile at position x (#Coils,Ny) 
    w = harm*pinv(cmap_tmp);                           % w        -> SMASH coil weights (af,#Coils) 
    sig_tmp = w*squeeze(sig_hb(:,:,x));                % sig_tmp  -> hyprid signal at position x (af, Ky_red)
    sig_out_hb(:,x) = sig_tmp(:);                      % sig_out_hb ->  reconstructed hyprid signals at position x (Ky_red*af,1)
end

sig_out = fftshift(fft(fftshift(sig_out_hb,2),[],2),2);     % fft hyprid signal (Ky,Nx) -> (Ky,Kx)  
img_out = ifftshift(ifft(ifftshift(sig_out_hb,1),[],1),1);  % fft hyprid signal (Ky,Nx) -> (Ny,Nx)

