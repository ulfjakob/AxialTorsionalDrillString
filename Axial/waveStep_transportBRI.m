function [x,y] = waveStep_transportBRI(p,x0,u)


    % Define convenience parameters
    % Drill string
    dx = p.L/p.P;
    dt = dx/max(p.c_a)*.99 ;  % As per the CFL cond.
    % BRI
    dxl = 1/p.Pl;
%     omegaMAX = max(p.omegaMAX,u.obEst);
    dtl = dxl/(p.omegaMAX*p.N/(2*pi))*.99;  % As per the CFL cond.
    
    n = max(ceil(p.dt/dt),ceil(p.dt/dtl));
    dt  = p.dt/n;

    % Parse input
    h = x0(1:p.P);              % Strain
    v = x0(p.P+1:2*p.P);        % Velocity
    l = x0(2*p.P+1:2*p.P+p.Pl); % Depth of cut

    % Compute Riemann invariants
    Aalpha = v+p.c_a*h; % Axial
    Abeta  = v-p.c_a*h;

    %%
    t = 0;
    for i=1:n
        % Topside BC
        Aalpha0  = -Abeta(1)*p.ka_L + 2*p.V0;
        
        %% Bit rock interaction
        % Bit velocity
        vb = 1/2*(Abeta(p.P)+Aalpha(p.P));
        
        ob = p.omega0;
        
        if ob>p.omegaMAX
%             keyboard;
            1;
        end

        % Combined depth of cut
        d = p.N * max( l(end) , 0);
        
        %% g function nonlinearity
        wc = d*p.K_a;   % Cutting component
        
        ga = 1/p.Wf * ( Aalpha(p.P)*p.Ap*p.E/p.c_a - wc);
        ga = (vb>0.0001) + ((vb)<=0.0001) * min( max(ga,0), 1);
        
       
        wf = p.Wf*ga;    % Contact component
        
        wb = wc + wf + u.wb;
        
        % BH BC
        AbetaPp1    = Aalpha(p.P)   - 2*wb*p.c_a/(p.Ap*p.E);
        
        %% 1st order Upwind
        % Drill string
        % Augment Riemann invariants with boundary values
        Aalpha_pad = [Aalpha0;  Aalpha];
        Abeta_pad  = [Abeta;    AbetaPp1];
% 
        Aalpha = Aalpha - dt/dx *p.c_a *diff(Aalpha_pad) -dt*p.ka*(Abeta+Aalpha);
        Abeta  = Abeta  + dt/dx *p.c_a *diff(Abeta_pad)  -dt*p.ka*(Abeta+Aalpha);
        
        % BRI
        % BRI BC is l(0) = 0
        l_pad = [0; l];
        l = l - dt/dxl *max(ob,0)*p.N/(2*pi) *diff(l_pad) + dt*vb;

        %% Second order Upwind
%         in1 = 1:p.P;
%         % Augment Riemann invariants with boundary values
%         Aalpha_pad = [2*Aalpha0-Aalpha(1); Aalpha0;  Aalpha];
%         Abeta_pad  = [Abeta;    AbetaPp1; 2*AbetaPp1-Abeta(p.P)];
% %         
%         Aalpha = Aalpha + dt/dx *p.c_a *1/2*(3*Aalpha_pad(in1)-4*Aalpha_pad(in1+1)+Aalpha_pad(in1+2))   -dt*p.ka*(Abeta+Aalpha);
%         Abeta  = Abeta  + dt/dx *p.c_a *1/2*(3*Abeta_pad(in1+2)-4*Abeta_pad(in1+1)+Abeta_pad(in1))      -dt*p.ka*(Abeta+Aalpha);
        
        % BRI BC is l(0) = 0
%         in1 = 1:p.Pl;
%         l_pad = [0-l(1); 0; l];
%         l = l - dt/dxl *max(ob,0)*p.N/(2*pi) *1/2*(3*l_pad(in1+2)-4*l_pad(in1+1)+l_pad(in1)) + dt*vb;

        t = t + dt;
    end
    %%
    v = 1/2         *(Abeta+Aalpha);
    h = 1/(2*p.c_a) *(Aalpha-Abeta);
    
    % Parse states
    x = nan(2*p.P+p.Pl,1);
    x(1:p.P)                = h;        % Strain
    x(p.P+1:2*p.P)          = v;        % Velocity
    x(2*p.P+1:2*p.P+p.Pl)   = l;        % Angular Velocity
    
    % Parse output
    y(1) = d;
    y(2) = 1-ga;
    y(3) = wf;
end




















