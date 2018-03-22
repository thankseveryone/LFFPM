function S=FPMdisplaced2(S,pixel_i,center1,center2,r,H_Fresnel,Im)

%O和P得到几何像的复振幅后，在频域中用CTF去截取，然后用用Fresnel函数

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
CTF = zeros(pixel_i);
CTF(xyc) = 1;
CTF = circshift(CTF,[-dr, -dc]);

F1 = fl.*H_Fresnel;
U1=ifft2(ifftshift(F1));

U2=angle(U1);
U3=sqrt(Im).*exp(1j*U2);

F3=fftshift(fft2(U3));
F4 = F3.*CTF;
phyipie = F4./H_Fresnel;

phyipie=circshift(phyipie,[dr, dc]);

S(xyc)=phyipie(xyc); 

end
% fl=zeros(pixel_i);
% regionc=zeros(pixel_i);
% 
% for i=1:pixel_i
%     for j=1:pixel_i
%         if (i-center1)^2+(j-center2)^2<=r^2;
%             fl(i,j)=S(i,j);
%             regionc(i,j)=1;
%         end
%     end
% end
% 
% dr=fix(center1-pixel_i/2-0.5);
% dc=fix(center2-pixel_i/2-0.5);
% fl=circshift(fl,[-dr, -dc]);  %这个时候应该得到了生成时的Fout
% 
% %应该为全零，检测过确实为全零。
% 
% % CTF=circshift(regionc, [-dr, -dc]);
% 
% U1=fl.*H_Fresnel;
% U2=ifft2(ifftshift(U1));
% 
% % U3=U2./abs(U2); 不能这样写！！！！！
% U3=angle(U2);
% U4=sqrt(Im).*exp(1j*U3); 
% %这样才是正确的振幅和相位！！！！！
% %这个时候得到的是对 离焦面 复振幅 的估计
% 
% U5=fftshift(fft2(U4)); 
% 
% % U5=U5.*CTF;
% 
% phyipie=U5./H_Fresnel;
% 
% phyipie=circshift(phyipie,[dr, dc]);
% 
% xyc=find(regionc==1);
% 
% S(xyc)=phyipie(xyc); 
% 
% end




