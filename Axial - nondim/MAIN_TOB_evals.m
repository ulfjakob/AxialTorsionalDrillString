%%
clc
clear;
close all;

%% Set parameter sets
Omega_vec = logspace(-.3,.7,100);
% Omega_vec = logspace(-.3,.7,20);
tN_vec    = 1./Omega_vec;
Ka_vec    = 20*ones(numel(tN_vec),1);
eta_vec   = .7*ones(numel(tN_vec),1);
MGVec = gainMargPipe(tN_vec,eta_vec,Ka_vec,'k');

%% Stability Map Location
% stabilityMap_pipe(tN_c,eta_c);

%% Simulation analysis
p = parameters;
% p.drawPeriod = 1200;
p.drawPeriod = inf;
T0 = 0;
T1 = 60*2;
% T1 = 10*1;
p.t = T0:p.dt:T1;
p.Nt = numel(p.t);

%%
p.eta   = 0.7;

% Ka_vec = [40];
Ka_vec = [40 20 10];
% Ka_vec = 10;
for K_a = Ka_vec
    p.K_a = K_a;
tic

p.k = -log(abs(p.eta))/2;
p.v0 = 1;
if p.eta<0
    p.w0 = p.Wf + 5;
end
hFig(1) = figure(3);
clf; hold off;
hFig(2) = figure(4);
clf; hold off;

% Log inits
TbAvg_vec = nan(size(tN_vec));
TcAvg_vec = nan(size(tN_vec));
TfAvg_vec = nan(size(tN_vec));
WbAvg_vec = nan(size(tN_vec));
WcAvg_vec = nan(size(tN_vec));
WfAvg_vec = nan(size(tN_vec));

VbAvg_vec = nan(size(tN_vec));


for i=1:numel(tN_vec)
    p.tN = tN_vec(i);
    p.omega = 1/(p.tN);
    
    [y,t] = runSim(p);
    
    %% Plottings
    ga = 1-y(2,1:p.Nt);
    d = y(1,1:p.Nt)*p.tN;
    tb = d + p.beta*(p.Wf*(ga-1));
    tc = d;                         % Cutting Torque
    tf = p.beta*(p.Wf*(ga-1));      % Contact Torque

    wc = d*p.K_a/p.ca;              % Cutting component
    wf = p.Wf*(ga-1);               % Contact component
    wb = wc + wf;
    
    vb = y(1,1:p.Nt);              % Bit velocity

    %%
%     sN1 = 1000;
%     sN2 = 5000;
    sN1 = 2000;
    sN2 = 11000;
    
    Tb_smooth = smooth(tb,sN1);
    TbAvg_vec(i) = Tb_smooth(sN2);
    Tc_smooth = smooth(tc,sN1);
    TcAvg_vec(i) = Tc_smooth(sN2);
    Tf_smooth = smooth(tf,sN1);
    TfAvg_vec(i) = Tf_smooth(sN2);
    
    Wb_smooth = smooth(wb,sN1);
    WbAvg_vec(i) = Wb_smooth(sN2);
    Wc_smooth = smooth(wc,sN1);
    WcAvg_vec(i) = Wc_smooth(sN2);
    Wf_smooth = smooth(wf,sN1);
    WfAvg_vec(i) = Wf_smooth(sN2);

    Vb_smooth = smooth(vb,sN1);
    VbAvg_vec(i) = Vb_smooth(sN2);

    %%
    if p.eta<0
        Vbar = 1./(p.k + p.K_a./Omega_vec)*(p.w0-p.Wf);
    else
        Vbar = 1;
    end
    Dbar = Vbar./Omega_vec;
    
    semilogx(hFig(1).Children,Omega_vec,TbAvg_vec,Omega_vec,TcAvg_vec,...
        Omega_vec,TfAvg_vec,Omega_vec,Dbar,'k--')
    semilogx(hFig(2).Children,Omega_vec,WbAvg_vec,Omega_vec,WcAvg_vec,...
    Omega_vec,WfAvg_vec,Omega_vec,VbAvg_vec)

    drawnow;
%     keyboard;
    toc
    i/numel(tN_vec)
end
% keyboard
%%
if p.eta<0
    Vbar = 1./(p.k + p.K_a./Omega_vec)*(p.w0-p.Wf);
else
    Vbar = 1;
end
Dbar = Vbar./Omega_vec;
figure(hFig(1)); clf;
semilogx(Omega_vec,TbAvg_vec,Omega_vec,TcAvg_vec,Omega_vec,TfAvg_vec)
hold all;
semilogx(Omega_vec,Dbar,'--k')
l = legend('Average torque on bit: $T_b-\beta W_f$','Average cutting torque: $D$',...
    'Average wearflat torque: $\beta W_f (g(V_b)-1)$','Equilibrium depth of cut: $\bar{D}$',...
    'location','best');
set(l,'interpreter','latex')
l = xlabel('$\Omega_b$'); set(l,'interpreter','latex');
tN = sprintf('$K_a$=%1.2f, $\\eta_a$ = %1.2f ',p.K_a,p.eta);
l = title(tN); set(l,'interpreter','latex');
grid on;
%%
figure(hFig(2));clf;
semilogx(Omega_vec,WbAvg_vec,Omega_vec,WcAvg_vec,Omega_vec,WfAvg_vec,Omega_vec,VbAvg_vec)
hold all;
% semilogx(Omega_vec,Dbar,'--k')
l = legend('Average weight on bit: $W_b- W_f$','Average cutting force: $\frac{K_a}{c_a}D$',...
    'Average wearflat force: $W_f (g(V_b)-1)$','Average bit velocity $V_b$', ...
    'location','best');
set(l,'interpreter','latex')
l = xlabel('$\Omega_b$'); set(l,'interpreter','latex');
tN = sprintf('$K_a$=%1.2f, $\\eta_a$ = %1.2f ',p.K_a,p.eta);
l = title(tN); set(l,'interpreter','latex');
grid on;
%%
if p.eta < 0
    fN1 = sprintf('AvgTorque%1.0fNEG',p.K_a);
    fN2 = sprintf('AvgWeight%1.0fNEG',p.K_a);
else
    fN1 = sprintf('AvgTorque%1.0f',p.K_a);
    fN2 = sprintf('AvgWeight%1.0f',p.K_a);
end
figure(hFig(1)); setLatexLabels2('x');
pdfmatlabfrag2(hFig(1),fN1);
figure(hFig(2)); setLatexLabels2('x');
pdfmatlabfrag2(hFig(2),fN2);
%%
end
%% Save PDF of figure
% print('-dpdf','LCLC')









