function [x,y] = waveStep_delayBRI(p,x0,u)

persistent XbHist PbHist

    % Define convenience parameters
    dx = p.L/p.P;
    dt = dx/max(p.c_a,p.c_a)*.99;  % As per the CFL cond.
    n = ceil(p.dt/dt)*8;
    dt = p.dt/n;

if isempty(XbHist)  || isempty(PbHist)
    XbHist = -((dt:dt:5).'*p.V0);
    PbHist = -((dt:dt:5).'*p.omega0);
end

    % Parse input
    h = x0(1:p.P);          % Strain
    v = x0(p.P+1:2*p.P);    % Velocity
    f = x0(2*p.P+1:3*p.P);  % Strain
    o = x0(3*p.P+1:4*p.P);  % Velocity
    Xb = x0(4*p.P+1);       % Bit Axial position
    Pb = x0(4*p.P+2);       % Bit Angular position
    
    % Compute Riemann invariants
    Aalpha = v+p.c_a*h; % Axial
    Abeta  = v-p.c_a*h;
    Talpha = o+p.c_t*f; % Torsional
    Tbeta  = o-p.c_t*f;

    %%
    t = 0;
    for i=1:n
        % Topside BC
        Aalpha0  = -Abeta(1) + 2*p.V0;
        Talpha0  = -Tbeta(1) + 2*p.omega0;
        
        %% Bit rock interaction
        % Bit velocity
        vb = 1/2*(Abeta(p.P)+Aalpha(p.P));
        ob = 1/2*(Tbeta(p.P)+Talpha(p.P));
        % Bit strain
        hb = 1/(2*p.c_a) *(Aalpha(p.P)-Abeta(p.P));

        % Find delay td_N = ind*dt
        ind = find( PbHist(1)-PbHist>2*pi/p.N ,1);
        if isempty(ind)
            error('History too short')
        end
        
        % Depth of cut
        d = max( XbHist(1)-XbHist( round(ind) ) , 0);
        
        % Weight/force required to keep bit stationary
        Ws = p.Ap*p.E/p.c_a * Aalpha(p.P);
%         if vb>0 && Ws>p.Wf
        if Ws>p.Wf
            wb = d*p.K_a + p.Wf;
        elseif Ws>0
            wb = Ws;
        else
            wb = d*p.K_a;
        end

        Ts = p.Jp*p.G/p.c_t * Talpha(p.P);
%         if ob>0 && Ts>p.Tf*p.Wf
        if Ts>p.Tf*p.Wf
            tb = d*p.K_t + p.Tf*p.Wf;
        elseif Ws>0
            tb = Ts;
        else
            tb = d*p.K_t;
        end

        % BH BC
        AbetaPp1    = Aalpha(p.P)   - 2*wb*p.c_a/(p.Ap*p.E);
        TbetaPp1    = Talpha(p.P)   - 2*tb*p.c_t/(p.Jp*p.G);
        
        % Augment Riemann invariants with boundary values
        Aalpha_pad = [Aalpha0; Aalpha];
        Abeta_pad  = [Abeta; AbetaPp1];
        Talpha_pad = [Talpha0; Talpha];
        Tbeta_pad  = [Tbeta; TbetaPp1];
        
        % 1st order Upwind 
        Aalpha = Aalpha - dt/dx *p.c_a *diff(Aalpha_pad) -p.ka/p.rho*(Abeta+Aalpha);
        Abeta  = Abeta  + dt/dx *p.c_a *diff(Abeta_pad)  -p.ka/p.rho*(Abeta+Aalpha);
        Talpha = Talpha - dt/dx *p.c_t *diff(Talpha_pad) -p.kt/p.rho*(Tbeta+Talpha);
        Tbeta  = Tbeta  + dt/dx *p.c_t *diff(Tbeta_pad)  -p.kt/p.rho*(Tbeta+Talpha);
        
        % Solve boundary ODE
        vb = 1/2*(Abeta(p.P)+Aalpha(p.P));
        ob = 1/2*(Tbeta(p.P)+Talpha(p.P));
        Xb = Xb + dt*vb;
        Pb = Pb + dt*ob;
        
        XbHist = [Xb; XbHist(1:end-1)];
        PbHist = [Pb; PbHist(1:end-1)];
        t = t + dt;
    end
    %%
    v = 1/2         *(Abeta+Aalpha);
    h = 1/(2*p.c_a) *(Aalpha-Abeta);
    o = 1/2         *(Tbeta+Talpha);
    f = 1/(2*p.c_t) *(Talpha-Tbeta);
    
    % Parse states
    x(1:p.P)        = h;        % Strain
    x(p.P+1:2*p.P)  = v;        % Velocity
    x(2*p.P+1:3*p.P)= f;        % Angular Strain
    x(3*p.P+1:4*p.P)= o;        % Angular Velocity
    x(4*p.P+1)      = Xb;       % Bit axial position
    x(4*p.P+2)      = Pb;       % Bit Angular position
    
    % Parse output
    y = d;
end




















