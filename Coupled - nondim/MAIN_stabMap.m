%%
clc
clear;
% close all;

%% Simulation analysis
p = parameters;

%% Sim stuff
T0 = 0;
T1 = 60*1;
p.t = T0:p.dt:T1;
p.Nt = numel(p.t);
p.eta_a = -0.70;
p.ka = -log(abs(p.eta_a));
p.v0 = 1;

KaVec = [10 20 40];
% KaVec = 20;
p.beta  = .1;

%%
nX = 30;
nY = 30;
V0Vec = logspace(-1,1.5,nY);
W0Vec = logspace(.3,3,nY);
% V0Vec = logspace(-.5,1.0,nY);
% W0Vec = logspace(.6,2,nY);
omMat = nan(nY,nX);
V0Mat = nan(nY,nX);
W0Mat = nan(nY,nX);
V0AvgMat = nan(nY,nX);
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
ii = 1;
for K_a = KaVec
    p.K_a = K_a;
    UTR_Mat = nan(nY,nX);
    nN = nY*nX;
    m = 1;
    for iY=1:nY
        for iX=1:nX
            p.tN0   = 1/omMat(iY,iX);
            if p.eta_a<0
                p.w0 = p.Wf + W0Mat(iY,iX);
                p.v0 = 0;
                V0Mat(iY,iX) = 0;
            else
                p.v0    = V0Mat(iY,iX);
            end

            plottings = 0;
            p.drawPeriod = 600;
            try
                [y,t] = runSim(p,plottings);
            catch me
                if strcmp(me.identifier,'upwind:CFLviolation')
                    UTR_Mat(iY,iX) = 1;
                else
                    UTR_Mat(iY,iX) = nan;
%                     rethrow(me)
                end
                m = m+1;
                continue;
            end
            vb_avg = smooth(y(4,1:p.Nt),round(p.Nt/2));
            V0AvgMat(iY,iX) = vb_avg( round(p.Nt*.9) );

            %%
%             figure(12);
%             clf;
%             plot(t,y(4,1:p.Nt))
%             hold on;
%             plot(t,y(5,1:p.Nt))
%             l = title(sprintf('$W_0 = %.1f, \\Omega_0 = %.2f$',...
%                 W0Mat(iY,iX),omMat(iY,iX)),'interpreter','latex');
% %             set(l,'interpreter','latex')
%             xlim([t(end)-10 t(end)])
% %             xlim([0 t(end)])
%             drawnow
            %%
            omega0 = 1/p.tN0;
            UTR = (omega0-min(y(5,p.Nt-20*round(1/p.dt):p.Nt)))/omega0;
            UTR_Mat(iY,iX) = UTR;
            m/nN

            m = m+1;
        end
    end
    %%
    h = figure; clf;
    if p.eta_a < 0
        fn = sprintf('plotSaves/UTR_K_a%d_beta%d_NEG',p.K_a,p.beta*100);
        UTR_contour_W(omMat,W0Mat,UTR_Mat,p)
        pdfmatlabfrag2(h,fn);
        
        fn2 = sprintf('plotSaves/UTR_K_a%d_beta%d_NEG_Vb',p.K_a,p.beta*100);
        h2 = figure; clf;
	
        pdfmatlabfrag2(h2,fn2);
    else
        fn = sprintf('plotSaves/UTR_K_a%d_beta%d',p.K_a,p.beta*100);
        UTR_contour(omMat,V0Mat,UTR_Mat,p)
        pdfmatlabfrag2(h,fn);
    end
    tS = sprintf('\\textbf{%s)}',char(96+ii));
    l = title(tS); set(l,'interpreter','latex','fontsize',20,...
    'Units', 'normalized','Position', [0.01 1], 'HorizontalAlignment', 'left');
    %%
    save(fn,'UTR_Mat','V0AvgMat')
    ii = ii+1;
end

%%
%
% plot([.25*KtVec+.05; .35*KtVec+3.5],KtVec,'r-')

%% Alternative plot
% 
% % nX = 30;
% % nY = 30;
% % V0Vec = logspace(-1,1.5,nY);
% omMat = nan(nY,nX);
% for i=1:nY
% %     omMat(i,:) = linspace(.22*KtVec(i)+.05,.35*KtVec(i)+3.5,nX);
% %     omMat(i,:) = logspace(log10(.22*KtVec(i)+.05),...
% %         log10(.37*KtVec(i)+4.5),nX);
%     omMat(i,:) = logspace(log10(.2*V0Vec(i)+.03),...
%         log10(1.2*V0Vec(i)+2.5),nX);
%     V0Mat(i,:) = ones(1,nX)*V0Vec(i);
% end
% 
% %
% p.K_a = 40;
% load('UTR_K_a40');
% 
% UTR_Mat(isnan(UTR_Mat))=1;
% 
% h = figure(3);
% clf
% % UTR_contour(omMat,V0Mat,UTR_Mat,p)
% % hold on;
% omMat1 = omMat;
% omMat2 = omMat;
% omMat3 = omMat;
% V0Mat1 = V0Mat;
% V0Mat2 = V0Mat;
% V0Mat3 = V0Mat;
% for iY=1:nY
%     for iX=1:nX
%         if UTR_Mat(iY,iX) > .95
%             omMat2(iY,iX) = nan;
%             V0Mat2(iY,iX) = nan;
%             omMat3(iY,iX) = nan;
%             V0Mat3(iY,iX) = nan;
%         elseif UTR_Mat(iY,iX) > .1
%             omMat1(iY,iX) = nan;
%             V0Mat1(iY,iX) = nan;
%             omMat3(iY,iX) = nan;
%             V0Mat3(iY,iX) = nan;
%         elseif UTR_Mat(iY,iX) >= 0
%             omMat1(iY,iX) = nan;
%             V0Mat1(iY,iX) = nan;
%             omMat2(iY,iX) = nan;
%             V0Mat2(iY,iX) = nan;
%         else 
%             omMat1(iY,iX) = nan;
%             V0Mat1(iY,iX) = nan;
%             omMat2(iY,iX) = nan;
%             V0Mat2(iY,iX) = nan;
%             omMat3(iY,iX) = nan;
%             V0Mat3(iY,iX) = nan;
%         end
%     end
% end
% % loglog(omMat1,V0Mat1,'xr')
% % loglog(omMat2,V0Mat2,'ob')
% % loglog(omMat3,V0Mat3,'.k')
% %
% loglog(nan,nan,'xr');
% hold on;
% loglog(nan,nan,'ob')
% loglog(nan,nan,'.k')
% 
% loglog(omMat1,V0Mat1,'xr')
% loglog(omMat2,V0Mat2,'ob')
% loglog(omMat3,V0Mat3,'.k')
% axis tight
% 
% l = ylabel('$V_0$'); set(l,'interpreter','latex');
% l = xlabel('$\Omega_0=1/t_N$'); set(l,'interpreter','latex');
% 
% tString = sprintf('$K_a=%d, \\ \\eta_a=%.2f \\ \\eta_t=%.2f$',p.K_a,p.eta_a,...
%     p.eta_t);
% l = title(tString); set(l,'interpreter','latex');
% 
% l = legend('$M_T>0.95$','$0.95\geq M_T> 0.1$','$0.1\geq M_T$',...
%     'location','northwest');
% set(l,'interpreter','latex');
% setLatexLabels2

% fn = sprintf('M_T_K_a%d',p.K_a);
% pdfmatlabfrag2(h,fn);
%%
% h = gcf;
% fn = sprintf('UTR_K_a%d',p.K_a);
% pdfmatlabfrag2(h,fn);

%%
% nX = 20;
% nY = 20;
% V0Vec = logspace(-.7,1.5,nY);
% omMat = nan(nY,nX);
% for i=1:nY
% %     omMat(i,:) = linspace(.22*KtVec(i)+.05,.35*KtVec(i)+3.5,nX);
% %     omMat(i,:) = logspace(log10(.22*KtVec(i)+.05),...
% %         log10(.37*KtVec(i)+4.5),nX);
%     omMat(i,:) = logspace(log10(.3*V0Vec(i)+.3),...
%         log10(.6*V0Vec(i)+4.5),nX);
%     V0Mat(i,:) = ones(1,nX)*V0Vec(i);
% end
% 
% %%



















