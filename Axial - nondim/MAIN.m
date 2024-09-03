%%
clc
clear;
% close all;


%% Simulation analysis
p = parameters;

T0 = 0;
T1 = 60*3;
p.t = T0:p.dt:T1;
p.Nt = numel(p.t);
p.drawPeriod = 150;
%%
p.eta = 0.7;
p.K_a = 40;
p.k = -log(abs(p.eta))/2;
% p.omega = 1.719;
% p.omega = 1.61;
p.omega = .65;
% p.w0 = p.Wf +(p.k + p.K_a/1 )*p.v0;
p.v0 = 1;
p.ca = 1.6;

%%
[y,t] = runSim(p);


%%
figure(11); clf;
tf = y(2,:)*p.Wf*p.beta;
Tf_smooth = smooth(tf,1000);
TfAvg_vec = Tf_smooth(5000)
plot(p.t,tf,p.t,Tf_smooth)

figure
plot(y(1,:))
hold on
plot(smooth(y(1,:),1e3))