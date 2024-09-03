%%
clc
clear;
close all;

%% Set parameter sets
tN_vec    = 1./[1 2 4 1 2 4 1 2 4];
% tN_vec    = sqrt(2)./[1 2 4 1 2 4 1 2 4];
eta_vec   = -ones(size(tN_vec))*0.7;
% Ka_vec    = ones(numel(tN_vec),1);
%% Stability Map Location
% stabilityMap_pipe(tN_vec,eta_vec);

%% Negative Gain Margins
% MGVec = gainMargPipe(tN_vec,eta_vec,Ka_vec,'k')
% Ka_vec = 1./MGVec.*Ka_vec;
% Ka_vec = Ka_vec.*[3 3 3 10 10 10 30 30 30].';
Ka_vec = [10 10 10 20 20 20 40 40 40]
MGVec = gainMargPipe(tN_vec,eta_vec,Ka_vec,'k')

%% Simulation analysis
p = parameters;
p.drawPeriod = inf;
% p.drawPeriod = 100;
T0 = 0;
T1 = 60*2;
% T1 = 10*1;
p.t = T0:p.dt:T1;
p.Nt = numel(p.t);

hFig = figure(10);
set(hFig, 'Position', [-900 0 560*2 420*2])
clf;
%%
for i=1:numel(tN_vec)
    p.tN = tN_vec(i);
    p.omega = 1/(p.tN);
    p.eta = eta_vec(i);
    p.k = -log(abs(p.eta))/2;
    p.K_a = Ka_vec(i);
    p.v0 = 0;
    
    if p.eta<0
        v0 = 1;
%         p.w0 = p.Wf +(p.k + p.K_a*1 )*v0 /2;
        p.w0 = p.Wf + 5;
%         p.w0 = 0;
    end

    [y,t] = runSim(p);
    
    %% Plottings
    figure(10);
    subplot(3,3,i)
    plot(p.t,y(5,1:p.Nt))
    hold all
    plot(p.t,y(1,1:p.Nt)*p.tN)
    plot(p.t,-y(2,1:p.Nt)*p.Wf)

    wb = y(1,1:p.Nt)*p.tN*p.K_a + y(3,1:p.Nt) - p.Wf;

    if i == 1
        l = legend('$V_b(\bar{t})$','$D(\bar{t})$','$W_f\big(g(\bar{t})-1\big)$',...
            'location','best');
        set(l,'interpreter','latex');
    end
    if i == 8
        l = xlabel('Nondimensional time $\bar{t}$'); set(l,'interpreter','latex');
    end
%     tS = sprintf('\\textbf{%s)} $\\Omega_0$=%1.2f, $K_a$=%1.0f, $K_AM_{G,a}$=%1.1f',char(96+i),...
%         p.omega,Ka_vec(i),MGVec(i));
    tS = sprintf('\\textbf{%s)} $\\Omega_0$=%1.2f, $K_a$=%1.0f, $\\eta_a$=%1.1f',char(96+i),...
        p.omega,Ka_vec(i),p.eta);
%     tString = sprintf('$\\textbf{1} $t_N$=%1.2f',MGVec(i));
    l = title(tS); set(l,'interpreter','latex');
    xlim(T1-[2 0])
    ylim([-5 10])
%     ylim([-6 6])
    
    %%
%     figure(21); clf;
%     plot(p.t,y(5,1:p.Nt))
%     hold all
%     plot(p.t,y(1,1:p.Nt)*p.tN)
%     plot(p.t,-y(2,1:p.Nt)*p.Wf)
% %     plot(p.t,wb)
% %     plot(p.t,y(4,1:p.Nt))
%     xlim(T1-[2 0])
% %     ylim([0 12])
%     ylim([-6 6])
%     grid on
    %%
    drawnow;
%     keyboard;
end
keyboard

%%
set(hFig,'Units','inches',...
    'Position',[1 0 7*1.72 5*1.72],...
    'PaperPositionMode','auto','PaperUnits','inches','PaperSize',[7*1.72 5*1.72])
%% Save PDF of figure
% print('-dpdf','LCLC_NEG')
% print('-dpdf','LCLC')

%%
keyboard


k = 1;
for k=1:9
% figure(2+k)
subplot(3,3,k)
% tString = sprintf('$t_N$ = %1.2f, $K_a$=%1.2f, $M_GK_a$=%1.2f',...
%     tN_vec(k),Ka_vec(k),MGVec(k));
% title(tString)
% xlim([58 60])
ylim([0 6])
end
%%
figure(hFig);


%%

% h = figure(4);
% pdfmatlabfrag2(h,'LC12');
% h = figure(5);
% pdfmatlabfrag2(h,'LC13');
% 
% h = figure(6);
% pdfmatlabfrag2(h,'LC21');
% h = figure(7);
% pdfmatlabfrag2(h,'LC22');
% h = figure(8);
% pdfmatlabfrag2(h,'LC23');
% 
% h = figure(9);
% pdfmatlabfrag2(h,'LC31');
% h = figure(10);
% pdfmatlabfrag2(h,'LC32');
% h = figure(11);
% pdfmatlabfrag2(h,'LC33');
%%
figure(11);
TbLp = zeros(size(t));
for k=2:numel(t)
    TbLp(k) = TbLp(k-1)*.98 + y(4,k-1)*.02;
end

clf;
plot(t,y(5,:))
hold on
plot(t,y(1,:))
plot(t,y(2,:)*p.Wf)
plot(t,y(4,:))
plot(t,TbLp(k),'-k','linewidth',2)
plot(t,t*0+1,'--k')
l = legend('$v_b(t)$','$d(t)$','$-\tilde{W}_f(t)$','$\tilde{T}_b(t)$',...
    '$\tilde{T}_b(t) Lowpass filter$',...
    'location','best');
set(l,'interpreter','latex');
title('Axial limit cycle')
l = xlabel('Nondimensional time $\bar{t}$'); set(l,'interpreter','latex');
l = ylabel('$v_b$ (m/s)'); set(l,'interpreter','latex');
ylabel('');
% xlim([55 60])
