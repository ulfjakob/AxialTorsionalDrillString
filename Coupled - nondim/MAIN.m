%%
clc
% clear;
close all;

%% Simulation analysis
p = parameters;
p.beta  = .1*1;

% p.K_a = 20;
% omegaVec = [0.4 0.6 0.9 1 1.8 3.2 2.2 3.5 6.3];
% v0Vec = [.4 .4 .4 2 2 2 10 10 10];

p.K_a = 20;
omegaVec = [1 1.8 3.2 2.2 3.5 6.3];
v0Vec = [2 2 2 10 10 10];

% omegaVec = [0.4 0.6 0.9 1 1.8 3.2 2.2 3.5 6.3];
% w0Vec = p.Wf + [15 15 15 20 20 20 35 35 35];

% p.K_a = 20;
% omegaVec = [.6 1.8 5 .6 1.8 5];
% w0Vec = p.Wf + [1 1 1 1.6 1.6 1.6]*40;

% w0Vec = p.Wf + 40;
% omegaVec = 4;

p.eta_a = 0.7;

%%
% omegaVec = logspace(-0.4,0.0,3);
% p.w0 = p.Wf + 4; p.eta_a = -0.7;
% omegaVec = logspace(-.4,0.5,3);
% p.w0 = p.Wf + 40; p.eta_a = -0.7;
% omegaVec = logspace(0.3,0.8,3);
% p.w0 = p.Wf + 50; p.eta_a = -0.7;

%%
% p.K_a = 20;
% omegaVec = 1.5;
% p.w0 = p.Wf + 25;
% p.eta_a = 0.7;
% p.v0 = 0.00;

% for ii=7:numel(omegaVec)
for ii=1
omega0 = omegaVec(ii);
if p.eta_a>0
    p.v0 = v0Vec(ii);
else
    p.w0 = w0Vec(ii);
    p.v0 = 0;
end


p.tN0   = 1/omega0;
p.ka = -log(abs(p.eta_a));

%% Sim stuff
T0 = 0;
T1 = 120;
p.t = T0:p.dt:T1;
p.Nt = numel(p.t);

plottings = 1;
p.drawPeriod = 200;

try
    [y,t] = runSim(p,plottings);
catch me
    if strcmp(me.identifier,'upwind:CFLviolation')
        keyboard;
    else
        rethrow(me)
    end
%     keyboard;
end

%%
Vb_Avg = smooth(y(4,1:p.Nt),round(p.Nt/3));
Vb_Avg = Vb_Avg(round(p.Nt*9/10))

% figure(12);
h = figure;
clf;
plot(t,y(4,1:p.Nt))
hold on;
plot(t,y(5,1:p.Nt))
tString = sprintf('$K_a=%d, \\ \\Omega_0=%.2f, \\ V_0=%.2f$',p.K_a,omega0,...
    p.v0);
% l = title(tString); set(l,'interpreter','latex');
tS = sprintf('\\textbf{%s)}',char(96+ii));
l = title(tS); set(l,'interpreter','latex','fontsize',20,...
    'Units', 'normalized','Position', [0.01 1], 'HorizontalAlignment', 'left');

l = legend('$V_b$','$\Omega_b$'); set(l,'interpreter','latex');
l = xlabel('$\bar{t}$'); set(l,'interpreter','latex');
% l = ylabel('$-$'); set(l,'interpreter','latex');
boldify
ylim([0 25])
% axis tight
xlim([t(end)-5 t(end)])
drawnow
%%
if p.eta_a<0
    fn = sprintf('plotSaves/coupledTS_Ka_%.0f_Omega_%.0f_W_%.0f_beta_%.0f_NEG',p.K_a,omega0*100,p.w0-p.Wf,p.beta*100);
else
    fn = sprintf('plotSaves/coupledTS_Ka_%.0f_Omega_%.0f_V_%.0f_beta_%.0f',p.K_a,omega0*100,p.v0*100,p.beta*100);
end
% set(h,'Units','inches',...
%     'Position',[1 3 7 5],...
%     'PaperPositionMode','auto','PaperUnits','inches','PaperSize',[7 5])
set(h,'Units','inches',...
    'Position',[1 3 4 3],...
    'PaperPositionMode','auto','PaperUnits','inches','PaperSize',[4 3])
print('-dpdf',fn)
%%
% pdfmatlabfrag2(h,fn);
end

















