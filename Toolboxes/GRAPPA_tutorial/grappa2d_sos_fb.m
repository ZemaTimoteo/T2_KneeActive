function [recon,g] = grappa2d_sos_fb(sig,acs,afy,afx,R); 

%   This is a 2D version of GRAPPA for acceleration in 2D (af = afy*afx)
%   In this code NO points in the readout direction are employed. Thus for
%   3D data sets GRAPPA has to be performed for every read-out point in
%   image space
%   
%   [recon,sigrecon,ws,ws_img] = grappa2d_fb(sig,acs,af);
%
%   This is a teaching version of the GRAPPA reconstruction code.
%   GRAPPA Weights are determined in k-space - Image reconstruction is
%   performed in image space (convolution theorem).
%       
%   IN:         sig                reduced data set         (#coils, Ky./afy, Kx./afx)
%               acs                autocalibration lines    (#coils, #acs_ky, #acs_ky)    
%               afy                Acceleration factor in y      (integer)
%               afx                Acceleration factor in x      (integer)
%               R                  noise covariance matrix  (#coils,#coils)    
%                                  
%
%   OUT:        recon              SOS combined GRAPPA image   (Ny, Nx)        
%               g                  SOS GRAPPA gfactor          (Ny, Nx)
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
%
%   
%   1) This code is strictly for non-commercial applications. The code is protected by
%      multiple patents.
%   2) This code is strictly for research purposes, and should not be used in any
%      diagnostic setting.
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
    error('Error! The number of coils has to be the same for both inputs!')
    
end

if afy*afx>nc
    error('Error! The acceleration factor must not exceed number of coils!')
      
end

if nargin < 5;
        R = eye(nc);
        disp('No noise correlations has been passed ......')
        disp('Assuming perfect coils without correlations......')
        no_corr = 1;
               
else
        sz = size(R);
        no_corr = 0;
        if size(sz)~=2 | sz(1)~=sz(2) | sz(1)~=nc
            error('Error! The dimension of noise covariance matrix has to be (#coils,#coils) !')
            
        end   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  calculate the weights for a 4x5 kernel in k-space from ACS data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                


sblx = 4*afx;                                          
sbly = 4*afy;                        
ypos = floor((sbly-afy)/2);    
xpos = floor((sblx-afx)/2);
   

idx = zeros(sbly,sblx);
idx(1:afy:end,1:afx:end)=1;

nsrc = size(idx(idx(:)==1),1);

idx=idx(:)==1;
src = zeros(nc*nsrc,(nyacs-sbly+1)*(nxacs-sblx+1));          
trg = zeros(nc*afy*afx,(nyacs-sbly+1)*(nxacs-sblx+1));   

cnt = 0;  % This is a lazy counter. could be done much faster.

for x=1:nxacs-sblx+1,
    for y=1:nyacs-sbly+1,
        cnt=cnt+1;
        
        tmp1 = acs(:,y:y+sbly-1,x:x+sblx-1);
        tmp2 = tmp1(:,idx);
        src(:,cnt) = tmp2(:); 
            
           
                % target points 
        tmp3 = acs(:,y+ypos:y-1+ypos+afy,x+xpos:x-1+xpos+afx);       
        trg(:,cnt) = tmp3(:);                  
            
                                                                                  
    end
end



ws=trg*pinv(src);

    
                                                                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                                   
%                                                                    
%  Calculate grappa weights for reconstruction in image space                                                                  
%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ws_k = zeros(nc,nc,ny*afy,nx*afx);             % prepare matrix for weights in zeropadded k-space 

ws_kernel_tmp = zeros(nc,afy,afx,nc,sbly*sblx);
ws_tmp = reshape(ws,[nc,afy,afx,nc,nsrc]);                                %Reshape weight set 

ws_kernel_tmp(:,:,:,:,idx(:)) = ws_tmp;
ws_kernel_tmp = flipdim(flipdim(reshape(ws_kernel_tmp,[nc afy afx nc sbly sblx]),5),6);

ws_kernel = zeros(nc,nc,sbly+afy-1,sblx+afx-1);


    for y = 1:afy,
        for x = 1:afx,
            ws_kernel(:,:,y:sbly+y-1,x:sblx+x-1) = ws_kernel(:,:,y:sbly+y-1,x:sblx+x-1) + squeeze(ws_kernel_tmp(:,y,x,:,:,:));    
        end
    end

[nc nc kszy kszx] = size(ws_kernel);

ws_k(:,:,ceil((ny-kszy)/2)+1:ceil((ny+kszy)/2),ceil((nx-kszx)/2)+1:ceil((nx+kszx)/2)) = ws_kernel;  

ws_img = ny*nx*ifftshift(ifftshift(ifft(ifft(ifftshift(ifftshift(ws_k,3),4),[],3),[],4),3),4);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Apply weights to folded images (GRAPPA reconstruction in image space -> Convolution Theorem)
%                                                                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


sig_red = zeros(nc,ny*afy,nx*afx);

sig_red(:,1:afy:end,1:afx:end) = sig;                                             
                                                                       
img_red = ifftshift(ifftshift(ifft(ifft(ifftshift(ifftshift(sig_red,2),3),[],2),[],3),2),3);

for k = 1:nc, 
 recon(k,:,:) = sum(squeeze(afy*afx*ws_img(k,:,:,:)).*img_red,1);             % This is the uncombined GRAPPA reconstruction;
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
        
            for y = 1:ny*afy,
                for x = 1:nx*afx,
                    tmp = squeeze(recon(:,y,x));             % GRAPPA-reco is used for calculating coil-combining coefficients n 
                                                             % ACS data can also be
                                                             % used
                    %sos = sqrt(abs(tmp'*inv(R)*tmp));       % SOS reconstruction
                    sos = sqrt(abs(tmp'*eye(nc)*tmp));
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
