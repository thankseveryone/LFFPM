%%
%说明：本程序用于仿真离焦情况的FPM。完成于2016年6月8日。

%强烈建议分节运行，即用 Ctrl+回车键 运行

%每一节开头均有内容说明。

%%
%加载物体和系统

load('USAFnew.mat');
size0=size(USAFnew,1);
USAF0=USAFnew;

%仿真系统使用参数如下
meterperpoint=0.000001; %物体大小1440微米

NA=0.055;

format long;
f=0.032; %物镜焦距 3.2cm
wuju=0.0470; %alpha=1平面物距为0.0465, 像距即微透镜面到主镜的距离是0.1026
%alpha=1 这样是放大2.2倍，是数值孔径匹配的

xiangju=1/(1/f-1/wuju);
M=xiangju/wuju;

alpha=xiangju/0.10262069;


lamda=0.0000005;
R=0.005; %主镜半径0.5厘米
k0=2*pi/lamda;
kradius=k0*NA;%频谱圈半径，单位 rad/m

image_size=size0*meterperpoint*M; %单位m

kperpoint=2*pi/0.00216; %在微透镜阵列面上，总尺寸为2160微米

rk=kradius/(M*kperpoint);%频谱圈半径所含像素数

f_microlens=20*0.000012;

%%
%计算菲涅尔衍射的滤波器 H_Fresnel
z=0.10262069-xiangju;

H_Fresnel_1440=zeros(size0);
for r=1:size0
    for c=1:size0
        fx=(r-0.5*size0-0.5)/0.00216;
        fy=(c-0.5*size0-0.5)/0.00216;
        s=(fx^2+fy^2)*(lamda^2);
        %FUout(r,c)=temp(r,c)*exp(-j*pi*lamda*z*(((r-0.5*sizer)/meterperpoint)^2+((c-0.5*sizec)/meterperpoint)^2));
        if s<1
            H_Fresnel_1440(r,c)=exp(1j*(2*pi/lamda)*z*sqrt(1-s));
        %fx, fy的值是 点数乘以1/size
        end
    end
end

%%
N=24; %共产生24张图，零频图在后面单独产生
%最终每张图是多大
sizeimage=1440;   %这个是我用于重聚的，插值后点数
% sizeimage=400;  %实验用的
incidentangle=zeros(1,N+1);
positionangle=zeros(1,N+1);
k=zeros(1,N+1);
kx=zeros(1,N+1);
ky=zeros(1,N+1);
kx_n=zeros(1,N+1);
ky_n=zeros(1,N+1);
centery=zeros(1,N+1);
centerx=zeros(1,N+1);
imageorigin=zeros(144,144,N+1);%存原始图的三维矩阵

%%
% 建立 LED 照明系统
% LED_h = 5.3; %单位cm 
% LED_d = 0.45; %单位cm %孔径变化时可直接替换数据，保证二者单位相同即可

LED_h=8.5;
LED_d=0.4;

ay_n = 2; %点的位置y坐标（-2，-1，0，1，2）
ax_n = 2; %点的位置x坐标（-2，-1，0，1，2）
% ay=-2, ax=-2表示最左上角
% x, y轴向右下为正。

n=0; %计数用

for ay=-ay_n:ay_n;%y轴坐标
    for ax=-ax_n:ax_n;%x轴坐标 
        if ay==0 && ax==0 
            continue %跳过零频图
        end
        n=n+1;
        if ax<0 %入射光的方位角（LED平面上的定位） 
            positionangle(n)=atan(ay/ax); 
        else
        	positionangle(n)=pi+atan(ay/ax);
        end
        
        incidentangle(n)= abs(atan(sqrt((LED_d*ax)^2+(LED_d*ay)^2)/LED_h));  
        %入射光的入射角（射向样本的角度）
        
        k(n)=k0*sin(incidentangle(n));
        kx(n)=k(n)*cos(positionangle(n));%求kx ky
        ky(n)=k(n)*sin(positionangle(n));
        
        %在像面上，空间频率要除以M
        kx_n(n)=kx(n)/(M*kperpoint);%转换为像素数，便于直接在矩阵中用 
        ky_n(n)=ky(n)/(M*kperpoint);        
        %这里的kx, ky就是我的cosiney*k0, cosinex*k0
        %kx, ky指向右下为正。
        
        centerx(n)=(sizeimage+1)/2+kx_n(n); %确定频谱中心
        centery(n)=(sizeimage+1)/2+ky_n(n);
        
    end
end

centerx(N+1)=(sizeimage+1)/2;
centery(N+1)=(sizeimage+1)/2; 

%%
%拍照，结果存储于zhengguoan变量
%本节自动计数

zhengguoan=zeros(1440,1440,49);

    %首先，经过缩放，得到几何光学像
    U_Gauss0=imresize(USAF0,image_size/0.002160);
    %你这个插值会带来高频，这就意味着你在Fresnel衍射中，得到的不是真实情况。
    cut_center=fix(size(U_Gauss0,1)/2+0.5);
    U_Gauss0=U_Gauss0(cut_center-719:cut_center+720,cut_center-719:cut_center+720);
    %至此，1440点表示2160微米。

for n=1:25
    
%     amplitude=cos(incidentangle(n))^2; 
    phase=zeros(size0); 
    for r=1:1440
        for c=1:1440
            phase(r,c)=exp(-1i*meterperpoint*1.5*((r-720.5)*ky(n)+(c-720.5)*kx(n))/M);
            % 查出这个和频谱圈不匹配，有负号才是对应的！！！！！
        end
    end
    
    U_Gauss=U_Gauss0.*phase;
    
    Fin=fftshift(fft2(U_Gauss)); %每格是1/0.00216
    %这个是标准物面的频谱
    
    %截止频率是:
    fxmax=(NA/(lamda*M)); %高斯像面截止频率    
    
    Fout=zeros(size0);
    for r=1:size0
        for c=1:size0
            fr=(r-0.5*(size0+1))/0.00216;
            fc=(c-0.5*(size0+1))/0.00216;
            if fr^2+fc^2<(fxmax)^2
                Fout(r,c)=Fin(r,c);
                Fin(r,c)=10000;
                %直接用像面截止频率算出的CTF的一格代表差不多500 m-1
                %物的频谱的一格代表差不多2000 m-1
                %所以直接相乘有问题
                %所以用*M来修正
            end
        end
    end

%     UArray=ifft2(ifftshift(Fout));%像面！
    
%     UArray=Fresnel(UArray,0.10262069-xiangju,lamda,0.00216/1440);

    F_CCD=Fout.*H_Fresnel_1440;
    U_CCD=ifft2(ifftshift(F_CCD));

    %这里直接是个CCD
    zhengguoan(:,:,n)=abs(U_CCD).^2;

    n
    fix(clock)
    %clear UonDet
end

%%
%观察位移
% for n=[7,8,9,13:15,17,18,19]
for n=1:25
figure(1); imagesc(squeeze(abs(zhengguoan(:,:,n)))); axis image; colormap(hot)
pause(1)
end

%%
%下面是FPM迭代。

guess=squeeze(zhengguoan(:,:,N+1)); %用零频图作为初始估计

%也可以用 S=zeros(sizeimage);做初始估计，但是会收敛得慢一些。

S=fftshift(fft2(guess));


%%
%迭代
%想要迭代几次，就反复运行这一节几次

for round=1
for j=[1:12,25,13:24]
    Im=zhengguoan(:,:,j);
    S=FPMdisplaced2(S,sizeimage,centery(j),centerx(j),rk,H_Fresnel_1440,Im);
    image=ifft2(ifftshift(S));
    figure(2); imagesc(abs(image)); axis image; colormap(hot);
    j 
end
end

