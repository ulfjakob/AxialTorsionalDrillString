%%
clc
clear;
% close all;

%
load('UTR_K_a60');

%% Physical Params
p.K_a       = 60;
p.eta_a     = .7;
p.eta_t     = .7;
p.cBar      = 1.6;


%% 

nX = 30;
nY = 30;
V0Vec = logspace(-1,1.5,nY);
omMat = nan(nY,nX);
for i=1:nY
%     omMat(i,:) = linspace(.22*KtVec(i)+.05,.35*KtVec(i)+3.5,nX);
%     omMat(i,:) = logspace(log10(.22*KtVec(i)+.05),...
%         log10(.37*KtVec(i)+4.5),nX);
    omMat(i,:) = logspace(log10(.2*V0Vec(i)+.03),...
        log10(1.2*V0Vec(i)+2.5),nX);
    V0Mat(i,:) = ones(1,nX)*V0Vec(i);
end



h = figure(2);
% clf
UTR_contour(omMat,V0Mat,UTR_Mat,p)

%% Nyquist stabMap
axStabMapForComp(p);

%% Alternative plot
figure(3)
% hold on;
omMat1 = omMat;
omMat2 = omMat;
omMat3 = omMat;
V0Mat1 = V0Mat;
V0Mat2 = V0Mat;
V0Mat3 = V0Mat;
for iY=1:nY
    for iX=1:nX
        if UTR_Mat(iY,iX) > .95
            omMat2(iY,iX) = nan;
            V0Mat2(iY,iX) = nan;
            omMat3(iY,iX) = nan;
            V0Mat3(iY,iX) = nan;
        elseif UTR_Mat(iY,iX) > .1
            omMat1(iY,iX) = nan;
            V0Mat1(iY,iX) = nan;
            omMat3(iY,iX) = nan;
            V0Mat3(iY,iX) = nan;
        elseif UTR_Mat(iY,iX) >= 0
            omMat1(iY,iX) = nan;
            V0Mat1(iY,iX) = nan;
            omMat2(iY,iX) = nan;
            V0Mat2(iY,iX) = nan;
        else 
            omMat1(iY,iX) = nan;
            V0Mat1(iY,iX) = nan;
            omMat2(iY,iX) = nan;
            V0Mat2(iY,iX) = nan;
            omMat3(iY,iX) = nan;
            V0Mat3(iY,iX) = nan;
        end
    end
end
% loglog(omMat1,V0Mat1,'xr')
% loglog(omMat2,V0Mat2,'ob')
% loglog(omMat3,V0Mat3,'.k')
%
loglog(omMat1,V0Mat1,'xr')
hold on;
loglog(omMat2,V0Mat2,'ob')
loglog(omMat3,V0Mat3,'.k')
axis tight

%%
% h = figure(2);
% pdfmatlabfrag2(h,'AxStabStab');














