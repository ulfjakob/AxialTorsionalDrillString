function [x,y] = waveStepUpwind(p,x0,u)


    % Define convenience parameters
    % Drill string
    dx = 1/p.P;
    dt = dx/max(1,p.ca)*.99 ;  % As per the CFL cond.
    % BRI
    dxl = 1/p.Pl;
    omegaMAX = 10/p.tN0;
    dtl = dxl/omegaMAX;  % As per the CFL cond.
    
    n = max(ceil(p.dt/dt),ceil(p.dt/dtl));
    dt  = p.dt/n;

    % Parse input
    h = x0(1:p.P);              % Strain
    v = x0(p.P+1:2*p.P);        % Velocity
    f = x0(2*p.P+1:3*p.P);      % Angular Strain
    o = x0(3*p.P+1:4*p.P);      % Angular Velocity
    l = x0(4*p.P+1:4*p.P+p.Pl); % Depth of cut
    
    % Compute Riemann invariants
    Aalpha = v+h*p.ca; % Axial
    Abeta  = v-h*p.ca;
    Talpha = o+f; % Torsional
    Tbeta  = o-f;

    %%
    t = 0;
    for i=1:n
        % Topside BC
        omega0 = 1/p.tN0;
        if isfield(p,'w0')
            Aalpha0  = Abeta(1) + 2*p.w0*p.ca;
        else
            Aalpha0  = -Abeta(1) + 2*p.v0;
        end
        Talpha0  = -Tbeta(1) + 2*omega0;
        
        %% Bit rock interaction
        % Bit velocity
        vb = 1/2*(Abeta(p.P)+Aalpha(p.P));
%         vb = p.v0;

        ob = 1/2*(Tbeta(p.P)+Talpha(p.P));
%         ob = 1/p.tN0;

        if ob>omegaMAX
            error('upwind:CFLviolation','ob larger than omegaMAX.');
        end

        % Normalized depth of cut
        d = max( l(end) , 0);
        
        %% g function nonlinearity
        wc = d*p.K_a;   % Axial cutting component
        tc = d;         % Torsional cutting component
        
        % Contact component
        % hb = 1/(2*p.ca) *(Aalpha(p.P)-Abeta(p.P));
        
        delta = .01;

        ga = 1/p.Wf * ( 1/(2*p.ca)*2*Aalpha(p.P) - wc - u.wb);
        ga = (vb>delta) + ((vb)<=delta) * min( max(ga,0), 1);

        wf = p.Wf*ga;
        tf = p.Wf*p.beta*ga;

        aStick = (ga<1) && (ga>0);
%         aStick = (vb<0.02)*(wres<0);
        
        % Torsional stick
        fb = 1/2 *(Talpha(p.P)-Tbeta(p.P));
        g_tt = fb - tc - p.Wf*p.beta;
        tStick = (ob<0.02) && (g_tt<0);

        tb = tc + tf + u.tb; % Total ToB
        wb = wc + wf + u.wb; % Total WoB

        % BH BC
        if tStick   % Zero velocity BC
            AbetaPp1    = -Aalpha(p.P);
            TbetaPp1    = -Talpha(p.P);
        elseif aStick
%             AbetaPp1    =  Aalpha(p.P)   - 2*wb*p.ca;
            AbetaPp1    = -Aalpha(p.P);
            TbetaPp1    =  Talpha(p.P)   - 2*tb;
        else
            AbetaPp1    =  Aalpha(p.P)   - 2*wb*p.ca;
            TbetaPp1    =  Talpha(p.P)   - 2*tb;
        end

        % Recompute bit velocity based on update stiction forces
%         vb = 1/2*(Abeta(p.P)+Aalpha(p.P));
%         ob = 1/2*(Tbeta(p.P)+Talpha(p.P));
        
        %% 1st order Upwind
        % Drill string
        % Augment Riemann invariants with boundary values
        Aalpha_pad = [Aalpha0;  Aalpha];
        Abeta_pad  = [Abeta;    AbetaPp1];
        Talpha_pad = [Talpha0;  Talpha];
        Tbeta_pad  = [Tbeta;    TbetaPp1];

        Aalpha = Aalpha - dt/dx *p.ca   *diff(Aalpha_pad) -dt*p.ka*(Abeta+Aalpha)/2;
        Abeta  = Abeta  + dt/dx *p.ca   *diff(Abeta_pad)  -dt*p.ka*(Abeta+Aalpha)/2;
        Talpha = Talpha - dt/dx         *diff(Talpha_pad) -dt*p.kt*(Tbeta+Talpha)/2;
        Tbeta  = Tbeta  + dt/dx         *diff(Tbeta_pad)  -dt*p.kt*(Tbeta+Talpha)/2;

        % Hack to ensure numeric stability of the inclusion
        if aStick || tStick
            Abeta(p.P)   = -Aalpha(p.P);
        end
        if tStick
            Tbeta(p.P) = - Talpha(p.P);
        end

        % BRI
        % BRI BC is l(0) = 0
        l_pad = [0; l];
        l = l - dt/dxl *max(ob,0) *diff(l_pad) + dt*vb;

        t = t + dt;
    end
    %%
    v = 1/2         *(Abeta+Aalpha);
    h = 1/(2*p.ca)  *(Aalpha-Abeta);
    o = 1/2         *(Tbeta+Talpha);
    f = 1/2         *(Talpha-Tbeta);
    
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




















