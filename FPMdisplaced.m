function O=FPMdisplaced(O,P,d_Fresnel,lamda,meterperpoint,Im)

%O是物体的复振幅透过率函数
u0=O.*P; 

fu0=fftshift(fft2(u0));
fu1=fu0.*CTF;

u1=ifft2(ifftshift(fu1));
%exit wave from the image plane

u2=Fresnel(u1,d_Fresnel,lamda,meterperpoint); %Propagate to the detector plane

u3=u2./abs(u2);

u4=u3.*sqrt(Im); %correct the modulus

u5=Fresnel(u4,-d_Fresnel,lamda,meterperpoint); %back-propagate to the image plane

delta_u=u5-u1; %两个exit wave之差

O=O+0.2*delta_u.*conj(P)/(max(max(abs(P)))^2); %更新透过率函数O


% dr=fix(center1-pixel_i/2-0.5);
% dc=fix(center2-pixel_j/2-0.5);
% fl=circshift(fl,[-dr, -dc]);
% 
% CTF=circshift(regionc, [-dr, -dc]);
% %FPM的核心
% 
% fl=fl.*H_Fresnel;  %这样得到C
% 
% Ui=ifft2(ifftshift(fl)); %这个时候得到的是对 离焦面 复振幅 的估计
% 
% 
% %精髓
% phase=Ui./abs(Ui);
% 
% Uiplus1=phase.*sqrt(Im);
% 
% phyipie=fftshift(fft2(Uiplus1)); %这就是他那里所说的C_hat
% 
% delta_G=(phyipie./H_Fresnel-fl./H_Fresnel).*abs(CTF).*conj(CTF)./(CTF.^2+0.01*ones(pixel_i));
% O=O+delta_G;
% 
% 
% phyipie=phyipie./H_Fresnel;
% phyipie=circshift(phyipie,[dr, dc]);
% 
% 
% xyc=find(regionc==1);
% O(xyc)=phyipie(xyc); 
