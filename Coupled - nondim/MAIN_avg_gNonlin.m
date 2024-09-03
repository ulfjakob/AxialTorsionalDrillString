%%
clc
clear;
% close all;

%% Simulation analysis
p = parameters;
ii = 1;

%%
fig_h1 = figure(21); clf;
set(fig_h1,'Units','inches','Position',[8 8 17 4.5]);

fig_h2 = figure(22); clf;
set(fig_h1,'Units','inches','Position',[7 7 17 4.5]);

%%
ii = 1;
for K_a = [10 20 40]
    tS = sprintf('plotSaves/UTR_K_a%d_beta10_NEG',K_a)

    load(tS);
    % load 'plotSaves/UTR_K_a10_beta10'
    % load 'plotSaves/UTR_K_a0'

    p.K_a = K_a;


    %% Sim stuff
    p.eta_a = -0.70;
    p.ka = -log(abs(p.eta_a));
    p.v0 = 1;

    % KaVec = [0 10 20 40];
    p.beta  = .1;

    %%
    nX = 30;
    nY = 30;
    V0Vec = logspace(-1,1.5,nY);
    W0Vec = logspace(.3,3,nY);
    omMat = nan(nY,nX);
    V0Mat = nan(nY,nX);
    W0Mat = nan(nY,nX);
    for i=1:nY
    %     omMat(i,:) = linspace(.22*KtVec(i)+.05,.35*KtVec(i)+3.5,nX);
    %     omMat(i,:) = logspace(log10(.22*KtVec(i)+.05),...
    %         log10(.37*KtVec(i)+4.5),nX);
        omMat(i,:) = logspace(log10(.2*V0Vec(i)+.03),...
            log10(1.2*V0Vec(i)+2.5),nX);
        V0Mat(i,:) = ones(1,nX)*V0Vec(i);
        W0Mat(i,:) = ones(1,nX)*W0Vec(i);
    end

    gBar = 1/p.Wf*( p.Wf + W0Mat - p.K_a*V0AvgMat./omMat);
    Wf = W0Mat - p.K_a*V0AvgMat./omMat;
    EvBar = V0AvgMat./(omMat.*W0Mat)*p.K_a;
    
    %%
    figure(fig_h1);
    subplot(1,3,ii)
    VbAvg_contour(omMat,W0Mat,EvBar,p)
    tS = sprintf('%s)',char(96+ii));
    l = title(tS); 
    set(l,'fontsize',20, ...
        'Units','normalized','Position', [0.01 1], 'HorizontalAlignment', 'left');
    %%
    figure(fig_h2);
    subplot(1,3,ii)
    UTR_contour_W(omMat,W0Mat,UTR_Mat,p)
    tS = sprintf('%s)',char(96+ii));
    l = title(tS); 
    set(l,'fontsize',20, ...
        'Units','normalized','Position', [0.01 1], 'HorizontalAlignment', 'left');
    ii = ii+1;
end

%%
figure(fig_h1);
gcf;
colorbar;
set(fig_h1,'Units','inches','Position',[8 8 17 4.5]);
ha = subplot(131);
hc = subplot(133);
set(hc,'position',[0.6916 0.1100 0.2134 0.8150]);


%%
filename = 'avgVb';
figure(fig_h1);
% hAxes=get(fig_h1,'CurrentAxes');
set(fig_h1,'Units','inches',...
    'Position',[8 8 17 4.5],...
    'PaperPositionMode','auto','PaperUnits','inches','PaperSize',[17 4.5])
print('-dpdf',filename)

filename = 'UTR_neg';
figure(fig_h2);
% hAxes=get(fig_h1,'CurrentAxes');
set(fig_h2,'Units','inches',...
    'Position',[8 8 17 4.5],...
    'PaperPositionMode','auto','PaperUnits','inches','PaperSize',[17 4.5])
print('-dpdf',filename)


































