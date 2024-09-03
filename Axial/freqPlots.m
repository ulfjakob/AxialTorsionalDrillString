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
Zca  =      1/p.c_a.*sqrt(1 + p.ka./(s));
Ga   = p.L *s/p.c_a.*sqrt(1 + p.ka./(s));
Zct  =      1/p.c_t.*sqrt(1 + p.kt./(s));
Gt   = p.L *s/p.c_t.*sqrt(1 + p.kt./(s));

ga_i = -1/p.c_a* 1./Zca.*tanh(Ga);
gt_i = -1/p.c_t* 1./Zct.*tanh(Gt);

% p.Mb = 10e3;
% % Collar and pipe system
% ga_Mb = p.Ap*sqrt(p.E*p.rho)*ga_i./( -ga_i*p.Mb.*s + p.Ap*sqrt(p.E*p.rho) );
% gt_Ib = p.Jp*sqrt(p.G*p.rho)*gt_i./( -gt_i*p.Ib.*s + p.Jp*sqrt(p.G*p.rho) );
%%
delay = (1-exp(-s*p.tN));

ga = ga_i;
gt = gt_i;

Ga =  -ga .*p.Ka./s.*delay;
Gt =   gt .*p.Kt./s.*delay;

%% BODE diagrams
pp = 1/2/pi;

% Subsystems
figure(1);clf
bodelin([ga gt],w);
l = legend('$g_a(s)$','$g_t(s)$');
set(l,'Interpreter','Latex');
subplot(211)
ylim([-40 40])

figure(2); clf
bodelin([Ga Gt],w); 
subplot(211); hold on;
plot(pp*w,zeros(size(w)),':k');
subplot(212); hold on;
plot(pp*w,-180*ones(size(w)),':k');
subplot(211)
l = legend('$G_a$','$G_t$','location','best');
set(l,'Interpreter','Latex');
title('Axial term: $G_a(s)$');
ylim([-40 50])

%%

h = figure(6); clf;
polarDiagram({Ga+Gt},h,30)
l = legend('$G$','location','best');
set(l,'Interpreter','Latex');




























