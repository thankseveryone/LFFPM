
function [S,fl,U3,phyipie]=FPMold(S,pixel_i,center1,center2,Im,r)

% S is the current best spectrum
% pixel_i is size(S,1)
% Im is the captured intensity of the current spectrum circle.
% r is the radius of the spectrum circle.

fl=zeros(pixel_i);

[rows, cols] = meshgrid(1:pixel_i);
distance_mat = sqrt((rows-center1).^2 + (cols-center2).^2);
xyc = find(distance_mat <=r);
fl(xyc) = S(xyc);

dr=fix(center2-pixel_i/2-0.5);
dc=fix(center1-pixel_i/2-0.5);
fl=circshift(fl,[-dr, -dc]);  %这个时候应该得到了生成时的Fout

U1=ifft2(ifftshift(fl));

U2=angle(U1);
U3=sqrt(Im).*exp(1j*U2);

phyipie=fftshift(fft2(U3));

phyipie=circshift(phyipie,[dr, dc]);

S(xyc)=phyipie(xyc); 



% imagesc(abs(S)); axis image; colormap(gray);
% system('pause');






           
    
