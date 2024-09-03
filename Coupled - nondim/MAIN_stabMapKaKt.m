%%
clc
clear;
% close all;

%% Simulation analysis
p = parameters;

%% Sim stuff
T0 = 0;
T1 = 60*3;
p.t = T0:p.dt:T1;
p.Nt = numel(p.t);
p.v0 = 0;

% KaVec = 20;
% KaVec = [10 15 25 40 60];
% KaVec = [15 25];
% p.K_a   = 20;
p.beta  = .5;
p.tN0   = .5;
%%
% nX = 14;
% nY = 7;
nX = 2;
nY = 2;
KtVec = logspace(0,1.3,nY);
KaMat = nan(nY,nX);
KtMat = nan(nY,nX);
for i=1:nY
    KaMat(i,:) = logspace(log10(KtVec(i)*.2),log10(KtVec(i)*10),nX);
    KtMat(i,:) = ones(1,nX)*KtVec(i);
end
%%
h = figure;
clf;
%%
UTR_Mat = nan(nY,nX);
nN = nY*nX;
m = 1;
for iY=1:nY
    for iX=1:nX
        p.K_a   = KaMat(iY,iX);
        p.K_t   = KtMat(iY,iX);

        plottings = 1;
        try
            [y,t] = runSim(p,plottings);
        catch me
            UTR_Mat(iY,iX) = nan;
            m = m+1;
            continue;
        end

        figure(12);
        clf;
        plot(t,y(4,1:p.Nt))
        hold on;
        plot(t,y(5,1:p.Nt))
%         xlim([t(end)-30 t(end)])
%         drawnow

        omega0 = 1/p.tN0;
        UTR = (omega0-min(y(5,p.Nt-20*round(1/p.dt):p.Nt)))/omega0;
        UTR_Mat(iY,iX) = UTR;
        m/nN

        m = m+1;
        
        %%
        figure(h); clf;
        UTR_Mat_c = min(UTR_Mat,1);
        contourf(KaMat,KtMat,UTR_Mat_c,[0 .1:.2:.9 1])
        hold on;
        plot(KaMat(iY,iX),KtMat(iY,iX),'xk')
        set(gca,'yscale','log')
        set(gca,'xscale','log')
        l = ylabel('$K_t$'); set(l,'interpreter','latex');
        l = xlabel('$K_a$'); set(l,'interpreter','latex');
        drawnow 
        
    end
end
%%
fn = sprintf('UTR_Om%d',1/p.tN0);
save(fn,'UTR_Mat')
%%


UTR_Mat_c = min(UTR_Mat,1);
contourf(KaMat,KtMat,UTR_Mat_c,[0 .1:.2:.9 1])
set(gca,'yscale','log')
set(gca,'xscale','log')
l = ylabel('$K_t$'); set(l,'interpreter','latex');
l = xlabel('$K_a$'); set(l,'interpreter','latex');

tString = sprintf('$\\omega=%d$',1/p.tN0);
l = title(tString); set(l,'interpreter','latex');

setLatexLabels2


%%
% pdfmatlabfrag2(h,fn);


%%
%
% plot([.25*KtVec+.05; .35*KtVec+3.5],KtVec,'r-')





















