% function freqPlots(p)

clc
% close all

p = parameters();


nw = 2e3;
% Angular frequency the Laplace variable is evaluated at
w = linspace(3e-3,10*2*pi,nw)';
% w = logspace(-1,2,nw)';
% Laplace variable
s = 1j*w;

%%
Zca  =    sqrt(1 + p.k./(s));
Ga   = s.*sqrt(1 + p.k./(s));

ga_i = -1./Zca.*tanh(Ga);
%%
delay = (1-exp(-s*p.tN));

ga = ga_i;
Ga =  -ga.*p.K_a./s.*delay;
%% BODE diagrams
pp = 1/2/pi;
opt.hz =1;

% Subsystems
figure(1);clf
bodelin([ga],w,opt);
l = legend('$g_a(s)$','$g_t(s)$');
set(l,'Interpreter','Latex');
subplot(211)
ylim([-40 40])

figure(2); clf
bodelin([Ga],w); 
subplot(211); hold on;
plot(pp*w,zeros(size(w)),':k');
subplot(212); hold on;
plot(pp*w,-180*ones(size(w)),':k');
subplot(211)
l = legend('$G_a$','location','best');
set(l,'Interpreter','Latex');
title('Axial term: $G_a(s)$');
ylim([-40 50])

%%

h = figure(6); clf;
polarDiagram({Ga},h,30)
l = legend('$G$','location','best');
set(l,'Interpreter','Latex');




























