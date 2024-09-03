function [x,y] = waveStepLFSG(p,x0,u)


% Define convenience parameters
% Drill string
dx = 1/p.P;
dt = dx/max(1)*.9 ;  % As per the CFL cond.
% BRI
dxl = 1/p.Pl;
%     omegaMAX = max(p.omegaMAX,u.obEst);
dtl = dxl/p.omegaMAX;  % As per the CFL cond.

n = max(ceil(p.dt/dt),ceil(p.dt/dtl));
dt  = p.dt/n;

% Parse input
h = x0(1:p.P);              % Strain
v = x0(p.P+1:2*p.P);        % Velocity
f = x0(2*p.P+1:3*p.P);      % Angular Strain
o = x0(3*p.P+1:4*p.P);      % Angular Velocity
l = x0(4*p.P+1:4*p.P+p.Pl); % Depth of cut

% Prealloc
h_new = nan(p.P,1);
v_new = nan(p.P,1);
f_new = nan(p.P,1);
o_new = nan(p.P,1);

%%
t = 0;
for i=1:n

    %% Bit rock interaction
    % Bit velocity
    vb = v(p.P);
    ob = o(p.P);
    ob = 1/p.tN0;

    if ob>p.omegaMAX
        keyboard;
        1;
    end

    % Normalized depth of cut
    d = max( l(end) , 0);

    %% g function nonlinearity
    wc = d*p.K_a;   % Axial cutting component
    tc = d*p.K_t;   % Torsional cutting component

    % Contact component
    ga = 1/p.Wf * ( h(p.P) - wc);
    ga = (vb>0.02) + ((vb)<=0.02) * min( max(ga,0), 1);
    aStick = (ga<1) * (ga>0);
    wf = p.Wf*ga;
    tf = p.Tf*ga;

    % Stick force
%         g_tt = Talpha(p.P)*p.tt - tc - tf;
%         g_tt = (ob<0.05) * min(g_tt,.00)  ;
%         if g_tt < 0
%            1;
%         end
    g_tt = 0;

%         g_wt = Aalpha(p.P) - wc - wf;
%         g_wt = g_wt * (g_tt<0);
    g_wt = 0;

    tb = tc + tf + g_tt + u.tb; % Total ToB
    wb = wc + wf + g_wt + u.wb; % Total WoB

    %% BC
    % Topside
    o(1) = 1/p.tN0;
    v(1) = p.v0;
    % Bottomhole
    h(p.P)    = wb;
    f(p.P)    = tb;

    %% Drill string: LFSG
    % Define sparse difference operators
    dtdx = dt/dx;
    UW  = gallery('tridiag',zeros(1,p.P-2),-dtdx*ones(1,p.P-1),dtdx*ones(1,p.P-2));

    UB = sparse(p.P-1,1,dtdx);
    UWu         = [UW,UB];

    in1 = 1:p.P-1;
    in2 = 2:p.P;
    if aStick
        UWo = gallery('tridiag',zeros(1,p.P-1),-dtdx*ones(1,p.P)  ,dtdx*ones(1,p.P-1));
        UWi = gallery('tridiag',zeros(1,p.P-3),-dtdx*ones(1,p.P-2),dtdx*ones(1,p.P-3));
        UBi = sparse(p.P-2,1,dtdx);
        UWBi = [UWi,UBi];
        ini = 2:p.P-1;
        
        v(p.P)  = 0;
        v(ini)  = v(ini)-       UWBi * h(in1) - dt*p.ka*v(ini);
        h       = h     -       UWo  * v;
    else
        v(in2) = v(in2) -       UWu * h - dt*p.ka*v(in2);
        h(in1) = h(in1) -       UWu * v;
    end
    o(in2) = o(in2) -1/p.tt*UWu * f - dt*p.ka*o(in2);
    f(in1) = f(in1) -       UWu * o;


    % BRI: Upwind
    % BRI BC is l(0) = 0
    l_pad = [0; l];
    l = l - dt/dxl *max(ob,0) *diff(l_pad) + dt*vb;

    t = t + dt;
end
%%

% Parse states
x(1:p.P)                = h;        % Strain
x(p.P+1:2*p.P)          = v;        % Velocity
x(2*p.P+1:3*p.P)        = f;        % Angular Strain
x(3*p.P+1:4*p.P)        = o;        % Angular Velocity
x(4*p.P+1:4*p.P+p.Pl)   = l;        % Angular Velocity

% Parse output
y(1) = d;
y(2) = 1-ga;
y(3) = wf;










