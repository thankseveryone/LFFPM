
function [regionc,cfimgm]=fnew(fms,pixel_i,pixel_j,center1,center2,Im);
r=29;
fl=zeros(pixel_i,pixel_j);
regionc=zeros(pixel_i,pixel_j);

for i=1:pixel_i;
    for j=1:pixel_j;  %��������Ƶ�׷�Χ
        if (i-center1)^2+(j-center2)^2<=r^2;
           fl(i,j)=fms(i,j);
           regionc(i,j)=1;  %�˲�λ��
        end
    end
end

%FPM�ĺ���

imgl=ifft2(ifftshift(fl));

imgml=imgl./abs(imgl); %���ȹ�һ��ֻ����λ

imgm=imgml.*Im; %Im�ǵ�ǰ�������ͼ�Ŀռ����ڿռ����Ƶ��Ȧ����λ���ü��ϡ���λ�õ�ǰ���½���ġ�

fimgm=fftshift(fft2(imgm));

cfimgm=fimgm.*regionc; %ֻ�������Ƶ��Ȧ�ڵ�ֵ





           
    
