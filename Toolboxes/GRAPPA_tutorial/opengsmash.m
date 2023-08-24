function [img_out, sig_out]=opengsmash(sig_in,cmap,af);

%
%   [img_out, sig_out]=opengsmash(sig_in,cmap,af);
%
%   G-SMASH reconstruction code
%       
%   IN:         sig_in             reduced k-space data     (#coils, Ky./af, Kx)
%               cmap               coil sensitivity maps    (#coils, Ny,Nx)
%               af                 Acceleration factor      (integer)
%
%   OUT:        img_out            Reconstructed images     (Ny, Nx)
%               sig_out            Reconstructed images     (ky, kx)
%
%   24.01.2006  Rita Nunes
%
%   modified by  Felix Breuer 04.06.2008 

[nc_cmap,ny,nx] = size(cmap);
[nc_sig_in,ky_red,kx] = size(sig_in);

if nc_cmap~=nc_sig_in
    disp('Error! The number of coils has to be the same for both inputs!')
    return;
else
    nc = nc_cmap;
end


if ky_red*af~=ny && ky_red~=ny || kx~=nx
    disp('Error! Size of folded data is not compatible with size of coil maps!')
    return;
end

cmap(isinf(cmap)) = 0;        %Just check to make sure there are no bad entries in the coil maps
cmap(isnan(cmap)) = 0;

tic;

hbd_sig = zeros(nc*ky_red,nx);
coef_mat = zeros(nc*ky_red,ny,nx);


sig_hbd  = ifftshift(ifft(ifftshift(sig_in,3),[],3),3);          % sig_hbd -> hybrid signal (#Coils,ky_red,Nx) -> (#Coils,ky_red,Nx)
cmap_hbd = fftshift(fft(fftshift(cmap,2),[],2),2);               % cmap_hbd -> hybrid sensitivities (#Coils,Ny,Nx) -> (#Coils,ky,Nx)

sig_hbd_mat = reshape(permute(sig_hbd,[2 1 3]),[ky_red*nc,nx]);  % sig_hbd_mat -> hybrid signal in 2D matrix

cmap_hbd_mat = zeros(ky_red*nc,ny,nx);                           % prepare hybrid sensitivity matrix

for m = 1:ky_red,
    cmap_hbd_mat(m:ky_red:end,:,:) = circshift(cmap_hbd,[0 ny/2-(m-1)*af 0]);    % cmap_hbd_mat ->  perform circular shift     
end

toc; 

tic;

for x=1:nx
   cmap_hbd_tmp = squeeze(cmap_hbd_mat(:,:,x));              % cmap_hbd_tmp -> hybrid sensitivity coefficients at position x   
   sig_hbd_tmp = squeeze(sig_hbd_mat(:,x));                  % sig_hbd_tmp  -> reduced hybrid data at position x
   sig_hbd_out(:,x)=pinv(cmap_hbd_tmp)*sig_hbd_tmp;          % sig_hbd_tmp  -> reconstructed hybrid data 
end

% (ky,x)->(y,x)
% we should be using ifft instead of fft, but then it would be necessary 
% to mirror the image along the y-axis and shift it by 1 pixel

img_out = ifftshift(fft(ifftshift(sig_hbd_out,1),[],1),1);

sig_out = fftshift(fft(fftshift(sig_hbd_out,2),[],2),2);
toc;



