
function [Pnew,Snew,Error]=fnew(S,P,pixel_i,pixel_j,center1,center2,Im,r)

center1=fix(center1);
center2=fix(center2);

move1=zeros(pixel_i, pixel_j);
move1(pixel_i-center1+1, pixel_j-center2+1)=1;

move2=zeros(pixel_i, pixel_j);
move2(center1, center2)=1;

Smove=conv2(S, move1, 'same');
phyi=Smove.*P;


%FPM的核心

Ui=ifft2(ifftshift(phyi));

Error=sum(sum(abs(abs(Ui)-abs(sqrt(Im)))))/sum(sum(abs(Ui)));

%imgml=Ui./(abs(Ui).^2+ones(pixel_i,pixel_j));
imgml=Ui./(abs(Ui)+0.1*ones(pixel_i,pixel_j));

%Uiplus1=imgml.*Im; %这是novel constraint
Uiplus1=imgml.*sqrt(Im);

phyipie=fftshift(fft2(Uiplus1));

% deltaS=conv2((phyipie-phyi).*conj(P)./(P.^2+0.01*ones(pixel_i,pixel_j)),move2,'same');
deltaS=conv2((phyipie-phyi).*conj(P)/(max(max(abs(P)))^2), move2,'same');

Snew=S+0.7*deltaS; %调整点三，参数alpha

% deltaP=(phyipie-phyi).*conj(Smove)./(Smove.^2+0.01*ones(pixel_i,pixel_j));
deltaP=(phyipie-phyi).*conj(Smove)/(max(max(abs(Smove)))^2);

Pnew=P+deltaP; %调整点四，参数beta



           
    
