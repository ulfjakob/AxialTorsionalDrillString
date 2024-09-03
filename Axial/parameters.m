%%
function p = parameters

% Physical parameters
p.ro = .0635;     %[m] Drill string inner radius
p.ri = .0540;     %[m] Drill string outer radius
p.rbo= .0762;     %[m] Drill collar outer radius
p.rbi= .0286;     %[m] Drill collar inner radius
p.a = .108;       %[m] Bit-rock interaction law parameter
p.epsilon = 60e6; %[Pa] Bit-rock interaction law parameter
p.zeta = .6;      %[-] Bit-rock interaction law parameter
p.sigma = 60e6;   %[Pa] Bit-rock contact stress
p.l = 11.e-3;     %[m] Total wear flat length
p.gamma = 1;      %[-] Bit geometry parameter
p.mu = .6;        %[-] Bit-rock friction coefficient
p.N = 4;          %[-] Number of blades of the PDC bit

% p.a = .108*4;       %[m] Bit-rock interaction law parameter
% p.epsilon = 6e6; %[Pa] Bit-rock interaction law parameter

p.Ap = pi*(p.ro^2-p.ri^2);    %[m^2] Drill string cross sectional area
p.Jp = pi/2*(p.ro^4-p.ri^4);  %[m^4] Drill string polar moment of inertia
p.Ab = pi*(p.rbo^2-p.rbi^2);  %[m^2] Drill collar cross sectional area
p.Jb = pi/2*(p.rbo^4-p.rbi^4);%[m^4] Drill collar polar moment of inertia

p.L = 2000;       %[m] Drill string length
p.E = 200e9;      %[Pa] Youngs modulus
p.G = 77e9;       %[m] Shear modulus
p.rho = 8e3;      %[kg/m3] Density

p.omega0 = 30/60*2*pi;   %[rad/s] Nominal angular velocity (found from RPM).
p.V0 = 30/3600;           %[m/s] Nominal rate of penetration (ROP)
p.dNom = (2*pi)*p.V0/(p.omega0*p.N);

p.ka   = 1 *.5;                      %[-] Axial damping
p.Mb = 12e3*0; %[kg] BHA mass

p.ka_L = 1;     % Reflection coefficient

%% computed quantities
p.c_a = sqrt(p.E/p.rho);      %[m/s] Axial wave velocity 
p.tN = 2*pi/(p.N*p.omega0);   %[s] Delay in the bit-rock interaction law

p.R    = p.V0/p.omega0;        %[-] ROP to RPM ration

%% BRI parameters
p.Ka = p.a*p.zeta*p.epsilon*p.N /(  p.Ap*sqrt(p.E*p.rho));
p.Kt =  p.R*p.a^2*p.epsilon*p.N /(2*p.Jp*sqrt(p.G*p.rho));

%% Misc params
p.P   = 2000;      % Number of cells in discretization
p.x   = linspace(0,p.L,p.P); 

% CFL considerations
p.dt = .001;    % Data sampling rate
dx = p.L/p.P;
dt = dx/max(p.c_a)*.99;  % As per the CFL cond.
n = ceil(p.dt/dt);
dt  = p.dt/n;

% BRI
p.omegaMAX = p.omega0*1; % for enforcing CFL condition
dxl = dt*(p.omegaMAX*p.N/(2*pi));
p.Pl = floor(1/dxl);
% p.Pl = 100;

p.Wf = p.a*p.l*p.sigma ;
p.Tf = p.gamma*p.mu*p.a;

p.K_a = p.a*p.zeta*p.epsilon *.1;
p.K_t = p.a^2*p.epsilon/2;

p.K_t = 0;
p.Tf = 0;
% p.K_a = 1e5*10*2;
% p.K_t = 9e5*.1*10;






















