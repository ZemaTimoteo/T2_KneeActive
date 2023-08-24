load cmap_head;

cmap = cmap_head;
img = phantom(128);
imagesc(abs(img));
nc=8

for k=1:nc,
    coil_img(k,:,:)=squeeze(cmap(k,:,:)).*img;
end


img_coil_noise=addnoise(coil_img,0.006);

imagesc(abs(squeeze(img_coil_noise(1,:,:))));
sig_coil_noise = fftshift(fftshift(fft(fft(fftshift(fftshift(img_coil_noise,2),3),[],2),[],3),2),3);
af=4;sx=16,sy=16
sig_red=sig_coil_noise(:,1:af:end,:);
imagesc(squeeze(abs(sig_coil_noise(1,:,:))))
acs = sig_coil_noise(:,65-sy:64+sy,65-sx:64+sx);
imagesc(squeeze(abs(acs(1,:,:))))
out = opengrappa(sig_red,acs,af,9,8);
subplot(2,2,1)
imagesc(squeeze(sum(abs(out).^2,1).^.5),[0 .4])
subplot(2,2,2)
out = opengrappa(sig_red,acs,af,7,6);
imagesc(squeeze(sum(abs(out).^2,1).^.5),[0 .4])
out = opengrappa(sig_red,acs,af,5,4);
subplot(2,2,3)
imagesc(squeeze(sum(abs(out).^2,1).^.5),[0 .4])
subplot(2,2,4)
out = opengrappa(sig_red,acs,af,3,2);
imagesc(squeeze(sum(abs(out).^2,1).^.5),[0 .4])
colormap gray