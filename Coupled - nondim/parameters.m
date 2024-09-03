function p = parameters

%% Parameters
p.ca = 1.6;

p.eta_a = 0.70;
p.ka = -log(p.eta_a);
p.eta_t = 0.70;
p.kt = -log(p.eta_t);

% p.P   = 500;      % Number of cells in discretization
p.P   = 1000;
p.x   = linspace(0,1,p.P);

% BC
p.tN0 = .25;
p.v0 = 1;

% CFL considerations
p.dt = .01;    % Data sampling rate
dx = 1/p.P;
dt = dx/max(1,p.ca)*.99 ;  % As per the CFL cond.
n = ceil(p.dt/dt);
dt  = p.dt/n;

% BRI
p.omegaMAX = 1/(p.tN0) * 2;
dxl = dt*p.omegaMAX;
p.Pl = floor(1/dxl);

p.Wf = 200;
p.beta = .1;
% p.Tf = 10;
% p.K_a = 20;
p.K_a = 20;

%% Sim stuff
T0 = 0;
T1 = 60*1;
p.t = T0:p.dt:T1;
p.Nt = numel(p.t);

