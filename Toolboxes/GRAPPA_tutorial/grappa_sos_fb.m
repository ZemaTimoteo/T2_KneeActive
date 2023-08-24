function [recon,g] = grappa_sos_fb(sig,acs,af,R); 

%  Felix Breuer 25/10/07
%  grappa weights transformed in image space 
%  grappa reconstruction takes place in image space
%
%  modified by Felix Breuer 10.05.08
%  grappa g-factor calculation included for sos-combined images
%   
%   [recon,sigrecon,ws,ws_img] = opengrappa(sig,acs,af);
%
%   This is a teaching version of the GRAPPA reconstruction code.
%   GRAPPA Weights are determined in k-space - Image reconstruction is
%   performed in image space
%       
%   IN:         sig                reduced data set         (#coils, Ky./af, Kx)
%               acs                autocalibration lines    (#coils, #acs_ky, #acs_ky)    
%               af                 Acceleration factor      (integer)
%               R                  noise covariance matrix  (#coils,#coils)    
%                                  
%
%   OUT:        recon              SOS combined GRAPPA image   (Ny, Nx)        
%               g                  GRAPPA gfactor              (Ny, Nx)
%               
%   Some things to think about when using this code:
%
%           -The ACS lines used for reconstruction are NOT included in the final reconstructed
%           data sets. If you want to do this, you can do it after reconstruction. Please check
%           that they are lined up with the reconstructed data. I have not checked this.
%
%           -Since the ACS lines are not included, feel free to use a different imaging sequence
%           for acquisition of the ACS lines. We have seen some advantage in doing this when
%           the sequence used for the reduced acquisition is very flow sensitive, for example.
%           This can also be faster in many cases.
%
%           -The 4x5 block size is normally optimal for most applications, but if you need  
%           something faster, you should try the 4x3 code which uses a smaller block size.
%
%   
%   1) This code is strictly for non-commercial applications. The code is protected by
%      multiple patents.
%   2) This code is strictly for research purposes, and should not be used in any
%      diagnostic setting.
%
%   22.10.2004  Mark Griswold (mark@physik.uni-wuerzburg.de)
%
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   This code also calculates the grappa g-factor for sos-combined GRAPPA images
%   directly from the GRAPPA reconstruction weights.
%   
%   The g-factor for each image pixel in each coil is simply:
%   g = sqrt(abs(diag(ws_img*R*ws_img')))./sqrt(diag(R))
%   
%   The noise covariance matrix R has to be derived from noise only samples (#coils,#samples)
%   R = cov(noise')
%
%   3)  The g-factor here is calculated for SOS combined images. Please note that the
%       g-factor for SOS is equal to the one for a normalized image
%       combination in the case of not too low SNR 
%   4)  It is very important to include potential noise correlations (R) for accurate
%       g-factor results!!!
%
%   10.05.2008 Felix Breuer (breuer@mr-bavaria.de)
%
%


[nc,ny,nx]=size(sig);
[nc_acs,nyacs,nxacs]=size(acs);     %Get the size of both the input data and the acs data

if nc_acs~=nc
    disp('Error! The number of coils has to be the same for both inputs!')
    recon_sos=[];
    g_sos =[];
    return
end

if af>nc
    disp('Error! The acceleration factor must not exceed number of coils!')
    recon_sos=[];
    g_sos =[];
    return    
end

if nargin < 4;
        R = eye(nc);
        disp('No noise correlations has been passed ......')
        disp('Assuming perfect coils without correlations......')
        no_corr = 1;
               
else
        sz = size(R);
        no_corr = 0;
        if size(sz)~=2 | sz(1)~=sz(2) | sz(1)~=nc
            disp('Error! The dimension of noise covariance matrix has to be (#coils,#coils) !')
            recon = [];
            g = [];
            return
        end   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  calculate the weights for a 4x5 kernel in k-space from ACS data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                

%   Simple example at 2x:            
%
%                   -O-O-O-O-O-     1   
%                   - - - - - -     2
%                   -O-O-O-O-O-     3       
%                   - - -X- - -     4         
%                   -O-O-O-O-O-     5
%                   - - - - - -     6
%                   -O-O-O-O-O-     7
%                   - - - - - -     8
%
%   The circles are the source points, and the X is the target point.


%  slide through the ACS data and collect all source points in matrix 
%  src and all target points in matrix trg 


srcx = 5;                         % should be odd                  
srcy = 4;                         % should be even 

src=zeros(nc*srcy*srcx,(nyacs-(srcy-1)*af)*(nxacs-(srcx-1)));
trg=zeros(nc*af,(nyacs-(srcy-1)*af)*(nxacs-(srcx-1)));


cnt = 0;  % This is a lazy counter. could be done much faster.

for xind=floor(srcx./2)+1:nxacs-floor(srcx./2),
    for yind=1:nyacs-(srcy-1)*af,
        cnt=cnt+1;
                
        src(:,cnt)=reshape(acs(:,yind:af:yind+(srcy-1)*af,xind-floor(srcx./2):xind+floor(srcx./2)),nc*srcy*srcx,1);                   
                                     
        
        targ(:,cnt)=reshape(acs(:,yind+(srcy./2-1)*af:yind+(srcy./2)*af-1,xind),nc*af,1);    %These are the target point
                                                                            
    end
end



ws=targ*pinv(src);

    
                                                                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                                   
%                                                                    
%  Calculate grappa weights for reconstruction in image space                                                                  
%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ws_k = zeros(nc,nc,ny*af,nx);

ws_tmp = reshape(ws,[nc,af,nc,srcy,srcx]);                                 %Reshape weight set 

ws_tmp = flipdim(flipdim(ws_tmp,4),5);                                     %flip source points in ky and kx for the convolution              


    for k=1:af,
        ws_kernel(:,:,k:af:af*srcy,:) = ws_tmp(:,k,:,:,:);                 %reconstruction kernel
    end    

ws_k(:,:,ceil((ny-srcy)*af/2)+1:ceil((ny+srcy)*af/2),ceil((nx-srcx)/2+1):ceil((nx+srcx)/2)) = ws_kernel;  %put reconstruction kernel in the center of matrix

ws_img = ny*nx*ifftshift(ifftshift(ifft(ifft(ifftshift(ifftshift(ws_k,3),4),[],3),[],4),3),4);    %Fouriertransform in image space and scale to final matrix size




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Apply weights to folded images (GRAPPA reconstruction in image space -> Convolution Theorem)
%                                                                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


sig_red = zeros(nc,ny*af,nx);

sig_red(:,1:af:end,:) = sig;                                             
                                                                       
img_red = ifftshift(ifftshift(ifft(ifft(ifftshift(ifftshift(sig_red,2),3),[],2),[],3),2),3);

for k = 1:nc, 
 recon(k,:,:) = sum(squeeze(af*ws_img(k,:,:,:)).*img_red,1);             % This is the uncombined GRAPPA reconstruction;
end

 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Calculate grappa g-factor for sos-combined images
%
%  See breuer et al. ISMRM 2008 abstract 10
%                                                                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargout == 2,
        if no_corr,
            disp('Be extremely careful! g-factor results might be inacurate......') 
        end
        
            disp('Calculating g-factor for SOS-combined GRAPPA images!')
            disp('GRAPPA Reconstruction is SOS-combined')
        
            for y = 1:ny*af,
                for x = 1:nx,
                    tmp = squeeze(recon(:,y,x));             % GRAPPA-reco is used for calculating coil-combining coefficients n 
                                                             % ACS data can also be
                                                             % used
                    sos = sqrt(abs(tmp'*tmp));        % SOS reconstruction
                   
                    n = tmp'./sos;                           % Coil combining coefficients
                    W = squeeze(ws_img(:,:,y,x));            % Weights in image space
                    recon_sos(y,x) = sos;
                    g(y,x) = sqrt(abs((n*W)*R*(n*W)'))./sqrt(abs((n*eye(nc))*R*(n*eye(nc))'));       % This is the generalized g-factor formulation 
                                                                                                         % for an arbitrary set of coil combining coefficients n
                end

            end
           recon = recon_sos;  
else
    disp('No g-factor is calculated!')
    disp('GRAPPA Reconstruction is uncombined')
end
