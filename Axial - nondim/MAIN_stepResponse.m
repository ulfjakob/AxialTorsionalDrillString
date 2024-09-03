%%
clc
clear;
% close all;

%% Simulation analysis
p = parameters;

p.eta = 0.70;
p.k = -log(p.eta)/2;
% eta = exp(k)*exp(-2)

p.K_a = 0;
p.Wf = 0;
p.v0 = 0;

T0 = 0;
T1 = 30;
p.t = T0:p.dt:T1;
p.Nt = numel(p.t);

[y,t] = runSim(p);

%%
t_a = 1;
figure(13);
clf;
plot(p.t,y(5,1:p.Nt)-1)
hold on;
% plot(p.t,1*exp(-p.k *p.t),'--k');
plot(p.t,1*exp(log(p.eta)/2 *p.t),'--k');
% plot(p.t,p.eta*exp(2)*exp(-p.t),'-.k');

% plot(p.t,p.t*0 +1/zeta_a,':k');
% plot([t_a t_a t_a*2 t_a*2 t_a*3 t_a*3 t_a*4 t_a*4 t_a*5 t_a*5 t_a*6 t_a*6]*2,...
%     [-1 1 1 -1 -1 1 1 -1 -1 1 1 -1],':k');
% plot(p.t,p.t*0 -1/zeta_a,':k');
plot(p.t,-exp(-p.k*p.t),'--k');
plot(p.t,p.t*0,'-k');
grid on;
l = legend('$\tilde{V}_b(t)$','$\pm e^{-k_at}$');
% l = legend('$V_b(t)$','$\pm \frac{e^{-k_at}}{\zeta_a}$',...1
%     '$2t_a,4t_a,6t_a,8t_a,\dots$');
set(l,'interpreter','latex','FontSize',12);
title('Bit velocity step response')
xlabel('Nondimensional time');
ylabel('$\tilde{V}_b$');
xlim([0 15])
% ylim([-8e-6 8e-6]);
%%
h = figure(13);
pdfmatlabfrag2(h,'stepResp');