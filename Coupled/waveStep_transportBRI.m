function [x,y] = waveStep_transportBRI(p,x0,u)


    % Define convenience parameters
    % Drill string
    dx = p.L/p.P;
    dt = dx/max(p.c_a,p.c_t)*.99;  % As per the CFL cond.
    % BRI
    dxl = 1/p.Pl;
%     omegaMAX = max(p.omegaMAX,u.obEst);
    dtl = dxl/p.omegaMAX;  % As per the CFL cond.
    
    n = max(ceil(p.dt/dt),ceil(p.dt/dtl));
    dt  = p.dt/n;

    % Parse input
    h = x0(1:p.P);              % Strain
    v = x0(p.P+1:2*p.P);        % Velocity
    f = x0(2*p.P+1:3*p.P);      % Strain
    o = x0(3*p.P+1:4*p.P);      % Velocity
    l = x0(4*p.P+1:4*p.P+p.Pl); % Depth of cut
    
    % Compute Riemann invariants
    Aalpha = v+p.c_a*h; % Axial
    Abeta  = v-p.c_a*h;
    Talpha = o+p.c_t*f; % Torsional
    Tbeta  = o-p.c_t*f;

    %%
    t = 0;
    for i=1:n
        % Topside BC
        Aalpha0  = -Abeta(1)*p.ka_L + 2*p.V0;
        Talpha0  = -Tbeta(1)*p.kt_L + 2*p.omega0;
        
        %% Bit rock interaction
        % Bit velocity
        vb = 1/2*(Abeta(p.P)+Aalpha(p.P));
        ob = 1/2*(Tbeta(p.P)+Talpha(p.P));
        
%         ob = p.omega0;
        
        if ob>p.omegaMAX
%             keyboard;
            1;
        end

        % Combined depth of cut
        d = p.N * max( l(end) , 0);
        
        %% g function nonlinearity
        wc = d*p.K_a;   % Cutting component
        tc = d*p.K_t;
        
        ga = 1/p.Wf * ( Aalpha(p.P)*p.Ap*p.E/p.c_a - wc);
        ga = (vb>0.001) + ((vb)<=0.001) * min( max(ga,0), 1);
        
        if ga == 0
            1;
        end
        
        wf = p.Wf*ga;    % Contact component
        tf = p.Tf*p.Wf*ga;
        
        % Stick force
        g_tt = Talpha(p.P)*p.Jp*p.G/p.c_t - tc - tf;
        g_tt = (ob<0.1) * min(g_tt,0);
        g_wt = Aalpha(p.P)*p.Ap*p.E/p.c_a - wc - wf;
        g_wt = g_wt * (abs(g_tt)>1e-3);
        
        tb = tc + tf + g_tt + u.tb;
        wb = wc + wf + g_wt + u.wb;

        % BH BC
        AbetaPp1    = Aalpha(p.P)   - 2*wb*p.c_a/(p.Ap*p.E);
        TbetaPp1    = Talpha(p.P)   - 2*tb*p.c_t/(p.Jp*p.G);
        
        %% 1st order Upwind
        % Drill string
        % Augment Riemann invariants with boundary values
        Aalpha_pad = [Aalpha0;  Aalpha];
        Abeta_pad  = [Abeta;    AbetaPp1];
        Talpha_pad = [Talpha0;  Talpha];
        Tbeta_pad  = [Tbeta;    TbetaPp1];

        Aalpha = Aalpha - dt/dx *p.c_a *diff(Aalpha_pad) -dt*p.ka*(Abeta+Aalpha);
        Abeta  = Abeta  + dt/dx *p.c_a *diff(Abeta_pad)  -dt*p.ka*(Abeta+Aalpha);
        Talpha = Talpha - dt/dx *p.c_t *diff(Talpha_pad) -dt*p.kt*(Tbeta+Talpha);
        Tbeta  = Tbeta  + dt/dx *p.c_t *diff(Tbeta_pad)  -dt*p.kt*(Tbeta+Talpha);
        
        % BRI
        % BRI BC is l(0) = 0
        l_pad = [0; l];
        l = l - dt/dxl *max(ob,0)*p.N/(2*pi) *diff(l_pad) + dt*vb;

        %% Second order Upwind
%         in1 = 1:p.P;
%         % Augment Riemann invariants with boundary values
%         Aalpha_pad = [2*Aalpha0-Aalpha(1); Aalpha0;  Aalpha];
%         Abeta_pad  = [Abeta;    AbetaPp1; 2*AbetaPp1-Abeta(p.P)];
%         Talpha_pad = [2*Talpha0-Talpha(1); Talpha0;  Talpha];
%         Tbeta_pad  = [Tbeta;    TbetaPp1; 2*TbetaPp1-Tbeta(p.P);];
%         
%         Aalpha = Aalpha + dt/dx *p.c_a *1/2*(3*Aalpha_pad(in1)-4*Aalpha_pad(in1+1)+Aalpha_pad(in1+2))   -dt*p.ka*(Abeta+Aalpha);
%         Abeta  = Abeta  + dt/dx *p.c_a *1/2*(3*Abeta_pad(in1+2)-4*Abeta_pad(in1+1)+Abeta_pad(in1))      -dt*p.ka*(Abeta+Aalpha);
%         Talpha = Talpha + dt/dx *p.c_t *1/2*(3*Talpha_pad(in1)-4*Talpha_pad(in1+1)+Talpha_pad(in1+2)) -dt*p.kt*(Tbeta+Talpha);
%         Tbeta  = Tbeta  + dt/dx *p.c_t *1/2*(3*Tbeta_pad(in1+2)-4*Tbeta_pad(in1+1)+Tbeta_pad(in1))  -dt*p.kt*(Tbeta+Talpha);
        
        % BRI BC is l(0) = 0
%         in1 = 1:p.Pl;
%         l_pad = [0-l(1); 0; l];
%         l = l - dt/dxl *max(ob,0)*p.N/(2*pi) *1/2*(3*l_pad(in1+2)-4*l_pad(in1+1)+l_pad(in1)) + dt*vb;

        t = t + dt;
    end
    %%
    v = 1/2         *(Abeta+Aalpha);
    h = 1/(2*p.c_a) *(Aalpha-Abeta);
    o = 1/2         *(Tbeta+Talpha);
    f = 1/(2*p.c_t) *(Talpha-Tbeta);
    
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
end




















