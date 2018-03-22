%%
%˵�������������ڷ����뽹�����FPM�������2016��6��8�ա�

%ǿ�ҽ���ֽ����У����� Ctrl+�س��� ����

%ÿһ�ڿ�ͷ��������˵����

%%
%���������ϵͳ

load('USAFnew.mat');
size0=size(USAFnew,1);
USAF0=USAFnew;

%����ϵͳʹ�ò�������
meterperpoint=0.000001; %�����С1440΢��

NA=0.055;

format long;
f=0.032; %�ﾵ���� 3.2cm
wuju=0.0470; %alpha=1ƽ�����Ϊ0.0465, ��༴΢͸���浽�����ľ�����0.1026
%alpha=1 �����ǷŴ�2.2��������ֵ�׾�ƥ���

xiangju=1/(1/f-1/wuju);
M=xiangju/wuju;

alpha=xiangju/0.10262069;


lamda=0.0000005;
R=0.005; %�����뾶0.5����
k0=2*pi/lamda;
kradius=k0*NA;%Ƶ��Ȧ�뾶����λ rad/m

image_size=size0*meterperpoint*M; %��λm

kperpoint=2*pi/0.00216; %��΢͸���������ϣ��ܳߴ�Ϊ2160΢��

rk=kradius/(M*kperpoint);%Ƶ��Ȧ�뾶����������

f_microlens=20*0.000012;

%%
%���������������˲��� H_Fresnel
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
        %fx, fy��ֵ�� ��������1/size
        end
    end
end

%%
N=24; %������24��ͼ����Ƶͼ�ں��浥������
%����ÿ��ͼ�Ƕ��
sizeimage=1440;   %������������ؾ۵ģ���ֵ�����
% sizeimage=400;  %ʵ���õ�
incidentangle=zeros(1,N+1);
positionangle=zeros(1,N+1);
k=zeros(1,N+1);
kx=zeros(1,N+1);
ky=zeros(1,N+1);
kx_n=zeros(1,N+1);
ky_n=zeros(1,N+1);
centery=zeros(1,N+1);
centerx=zeros(1,N+1);
imageorigin=zeros(144,144,N+1);%��ԭʼͼ����ά����

%%
% ���� LED ����ϵͳ
% LED_h = 5.3; %��λcm 
% LED_d = 0.45; %��λcm %�׾��仯ʱ��ֱ���滻���ݣ���֤���ߵ�λ��ͬ����

LED_h=8.5;
LED_d=0.4;

ay_n = 2; %���λ��y���꣨-2��-1��0��1��2��
ax_n = 2; %���λ��x���꣨-2��-1��0��1��2��
% ay=-2, ax=-2��ʾ�����Ͻ�
% x, y��������Ϊ����

n=0; %������

for ay=-ay_n:ay_n;%y������
    for ax=-ax_n:ax_n;%x������ 
        if ay==0 && ax==0 
            continue %������Ƶͼ
        end
        n=n+1;
        if ax<0 %�����ķ�λ�ǣ�LEDƽ���ϵĶ�λ�� 
            positionangle(n)=atan(ay/ax); 
        else
        	positionangle(n)=pi+atan(ay/ax);
        end
        
        incidentangle(n)= abs(atan(sqrt((LED_d*ax)^2+(LED_d*ay)^2)/LED_h));  
        %����������ǣ����������ĽǶȣ�
        
        k(n)=k0*sin(incidentangle(n));
        kx(n)=k(n)*cos(positionangle(n));%��kx ky
        ky(n)=k(n)*sin(positionangle(n));
        
        %�������ϣ��ռ�Ƶ��Ҫ����M
        kx_n(n)=kx(n)/(M*kperpoint);%ת��Ϊ������������ֱ���ھ������� 
        ky_n(n)=ky(n)/(M*kperpoint);        
        %�����kx, ky�����ҵ�cosiney*k0, cosinex*k0
        %kx, kyָ������Ϊ����
        
        centerx(n)=(sizeimage+1)/2+kx_n(n); %ȷ��Ƶ������
        centery(n)=(sizeimage+1)/2+ky_n(n);
        
    end
end

centerx(N+1)=(sizeimage+1)/2;
centery(N+1)=(sizeimage+1)/2; 

%%
%���գ�����洢��zhengguoan����
%�����Զ�����

zhengguoan=zeros(1440,1440,49);

    %���ȣ��������ţ��õ����ι�ѧ��
    U_Gauss0=imresize(USAF0,image_size/0.002160);
    %�������ֵ�������Ƶ�������ζ������Fresnel�����У��õ��Ĳ�����ʵ�����
    cut_center=fix(size(U_Gauss0,1)/2+0.5);
    U_Gauss0=U_Gauss0(cut_center-719:cut_center+720,cut_center-719:cut_center+720);
    %���ˣ�1440���ʾ2160΢�ס�

for n=1:25
    
%     amplitude=cos(incidentangle(n))^2; 
    phase=zeros(size0); 
    for r=1:1440
        for c=1:1440
            phase(r,c)=exp(-1i*meterperpoint*1.5*((r-720.5)*ky(n)+(c-720.5)*kx(n))/M);
            % ��������Ƶ��Ȧ��ƥ�䣬�и��Ų��Ƕ�Ӧ�ģ���������
        end
    end
    
    U_Gauss=U_Gauss0.*phase;
    
    Fin=fftshift(fft2(U_Gauss)); %ÿ����1/0.00216
    %����Ǳ�׼�����Ƶ��
    
    %��ֹƵ����:
    fxmax=(NA/(lamda*M)); %��˹�����ֹƵ��    
    
    Fout=zeros(size0);
    for r=1:size0
        for c=1:size0
            fr=(r-0.5*(size0+1))/0.00216;
            fc=(c-0.5*(size0+1))/0.00216;
            if fr^2+fc^2<(fxmax)^2
                Fout(r,c)=Fin(r,c);
                Fin(r,c)=10000;
                %ֱ���������ֹƵ�������CTF��һ�������500 m-1
                %���Ƶ�׵�һ�������2000 m-1
                %����ֱ�����������
                %������*M������
            end
        end
    end

%     UArray=ifft2(ifftshift(Fout));%���棡
    
%     UArray=Fresnel(UArray,0.10262069-xiangju,lamda,0.00216/1440);

    F_CCD=Fout.*H_Fresnel_1440;
    U_CCD=ifft2(ifftshift(F_CCD));

    %����ֱ���Ǹ�CCD
    zhengguoan(:,:,n)=abs(U_CCD).^2;

    n
    fix(clock)
    %clear UonDet
end

%%
%�۲�λ��
% for n=[7,8,9,13:15,17,18,19]
for n=1:25
figure(1); imagesc(squeeze(abs(zhengguoan(:,:,n)))); axis image; colormap(hot)
pause(1)
end

%%
%������FPM������

guess=squeeze(zhengguoan(:,:,N+1)); %����Ƶͼ��Ϊ��ʼ����

%Ҳ������ S=zeros(sizeimage);����ʼ���ƣ����ǻ���������һЩ��

S=fftshift(fft2(guess));


%%
%����
%��Ҫ�������Σ��ͷ���������һ�ڼ���

for round=1
for j=[1:12,25,13:24]
    Im=zhengguoan(:,:,j);
    S=FPMdisplaced2(S,sizeimage,centery(j),centerx(j),rk,H_Fresnel_1440,Im);
    image=ifft2(ifftshift(S));
    figure(2); imagesc(abs(image)); axis image; colormap(hot);
    j 
end
end

