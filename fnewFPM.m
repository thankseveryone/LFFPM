
function [regionc,cfimgm]=fnew(fms,pixel_i,pixel_j,center1,center2,Im);
r=29;
fl=zeros(pixel_i,pixel_j);
regionc=zeros(pixel_i,pixel_j);

for i=1:pixel_i;
    for j=1:pixel_j;  %遍历整个频谱范围
        if (i-center1)^2+(j-center2)^2<=r^2;
           fl(i,j)=fms(i,j);
           regionc(i,j)=1;  %滤波位置
        end
    end
end

%FPM的核心

imgl=ifft2(ifftshift(fl));

imgml=imgl./abs(imgl); %幅度归一，只留相位

imgm=imgml.*Im; %Im是当前这幅低清图的空间域。在空间域把频谱圈的相位作用加上。相位用当前最新结果的。

fimgm=fftshift(fft2(imgm));

cfimgm=fimgm.*regionc; %只保留这个频谱圈内的值





           
    
