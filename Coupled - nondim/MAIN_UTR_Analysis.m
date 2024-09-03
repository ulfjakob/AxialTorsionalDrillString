%%
clc
clear;
% close all;

%% Simulation analysis
p = parameters;

load 'plotSaves/UTR_K_a10_NEG'
% load 'plotSaves/UTR_K_a10_beta10'
% load 'plotSaves/UTR_K_a0'

p.K_a = 10;
ii = 1;

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

%%
h = figure(20+ii); clf;
fn = sprintf('plotSaves/UTR_K_a%d_NEG',p.K_a);
UTR_contour_W(omMat,W0Mat,UTR_Mat,p)
tS = sprintf('\\textbf{%s)}',char(96+ii));
l = title(tS); set(l,'interpreter','latex','fontsize',20,...
    'Units', 'normalized','Position', [0.01 1], 'HorizontalAlignment', 'left');

% h = figure(30+ii); clf;
% UTR_contour(omMat,V0Mat,UTR_Mat,p)
% tS = sprintf('\\textbf{%s)}',char(96+ii));
%     l = title(tS); set(l,'interpreter','latex','fontsize',20,...
%     'Units', 'normalized','Position', [0.01 1], 'HorizontalAlignment', 'left');
% ylim([.1 10^1.5])
% UTR_contour(omMat,V0Mat,UTR_Mat,p)    
% hold on;
% plot([.1 10], 40*[1 1],'--r')

% hold on;
% plot(smooth(OmLine,1),V0Vect,'-r','linewidth',3)

%%
boldify

tS = sprintf('\\textbf{%s)}',char(96+ii));
l = title(tS); set(l,'interpreter','latex','fontsize',20,...
'Units', 'normalized','Position', [0.01 1], 'HorizontalAlignment', 'left');

%%
h = figure(23);
fn = sprintf('plotSaves/UTR_K_a%d_beta%d',p.K_a,p.beta*100);
% fn = sprintf('plotSaves/UTR_K_a%d_beta%d_NEG',p.K_a,p.beta*100);
% fn2 = sprintf('plotSaves/UTR_K_a%d_beta%d_NEG_Vb',p.K_a,p.beta*100);
% pdfmatlabfrag2(h,fn);





