function p = parameters

%% Parameters
p.tN = 0.6;
p.eta = 0.70;
p.k = -log(p.eta)/2;
p.ca = 1.6;
% Torque estimate param
p.beta = .1;

p.P   = 500;      % Number of cells in discretization
p.x   = linspace(0,1,p.P);

% BC
p.v0 = 1;

% CFL considerations
% p.dt = .05;    % Data sampling rate
p.dt = .01;    % Data sampling rate
dx = 1/p.P;
dt = dx/max(1,p.ca)*.99 ;  % As per the CFL cond.

n = ceil(p.dt/dt);
dt  = p.dt/n;

% BRI
p.omega = 1/(p.tN);
dxl = dt*p.omega;
p.Pl = floor(1/dxl);

p.Wf = 100;
p.K_a = 20;

% To change topside BC:
% p.w0 = p.Wf +(p.k + p.K_a*p.tN )*p.v0;

%% Sim stuff
T0 = 0;
T1 = 60*1;
p.t = T0:p.dt:T1;
p.Nt = numel(p.t);

p.drawPeriod = 5000;