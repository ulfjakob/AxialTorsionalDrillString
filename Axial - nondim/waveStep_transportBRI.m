function [x,y] = waveStep_transportBRI(p,x0,u)


    % Define convenience parameters
    % Drill string
    dx = 1/p.P;
    dt = dx/max(p.ca)*.99*1;  % As per the CFL cond.
    % BRI
    dxl = 1/p.Pl;
%     omegaMAX = max(p.omegaMAX,u.obEst);
    dtl = dxl/(p.omega)*.99;  % As per the CFL cond.
    
    n = max(ceil(p.dt/dt),ceil(p.dt/dtl));
    dt  = p.dt/n;

    % Parse input
    h = x0(1:p.P);              % Strain
    v = x0(p.P+1:2*p.P);        % Velocity
    l = x0(2*p.P+1:2*p.P+p.Pl); % Depth of cut

    % Compute Riemann invariants
    Aalpha = v+h*p.ca; % Axial
    Abeta  = v-h*p.ca;
    
    % Average Torque on bit
    TbCum = 0;

    %%
    t = 0;
    for i=1:n
        % Topside BC
        if isfield(p,'w0')
            Aalpha0  = Abeta(1) + 2*p.w0*p.ca;
        else
            Aalpha0  = -Abeta(1) + 2*p.v0;
        end
                
        %% Bit rock interaction
        % Bit velocity
        vb = 1/2*(Abeta(p.P)+Aalpha(p.P));
        
        % Combined depth of cut
        d = max( l(end) , 0);
        
        %% g function nonlinearity
        wc = d*p.K_a/p.ca;   % Cutting component
        
        delta = .01;
        
        ga = 1/p.Wf * ( 1/(2*p.ca)*2*Aalpha(p.P) - wc - u.wb);
        ga = (vb>delta) + ((vb)<=delta) * min( max(ga,0), 1);
%         ga = (vb>0.02) + ((vb)<=0.02) * min( max(ga,0), 1);
%         ga = (vb>0.05) + ((vb)<=0.05) * min( max(ga,0), 1);
%         ga = (vb>0.10) + ((vb)<=0.10) * min( max(ga,0), 1);
        
        wf = p.Wf*ga;           % Contact component
        wb = wc + wf + u.wb;    % Total WOB
        
        % Torque on bit
        tb = d + p.beta*(p.Wf*(ga-1));
        
        % BH BC
        AbetaPp1    = Aalpha(p.P)   - 2*wb*p.ca;
        
        %% 1st order Upwind
        % Drill string
        % Augment Riemann invariants with boundary values
        Aalpha_pad = [Aalpha0;  Aalpha];
        Abeta_pad  = [Abeta;    AbetaPp1];
% 
        Aalpha = Aalpha - dt/dx*p.ca *diff(Aalpha_pad) -dt*p.k*(Abeta+Aalpha);
        Abeta  = Abeta  + dt/dx*p.ca *diff(Abeta_pad)  -dt*p.k*(Abeta+Aalpha);
        
        % Hack to ensure numeric stability of the inclusion
        if ga<1 && ga>0
            Abeta(p.P)   = -Aalpha(p.P);
        end
        
        % BRI
        % BRI BC is l(0) = 0
        l_pad = [0; l];
        l = l - dt/dxl *p.omega *diff(l_pad) + dt*vb;

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
        
        %%
        % Log averaged torque on bit
        TbCum = TbCum + tb;
    end
    %%
    v = 1/2         *(Abeta+Aalpha);
    h = 1/(2*p.ca)  *(Aalpha-Abeta);
    
    % Parse states
    x = nan(2*p.P+p.Pl,1);
    x(1:p.P)                = h;        % Strain
    x(p.P+1:2*p.P)          = v;        % Velocity
    x(2*p.P+1:2*p.P+p.Pl)   = l;        % Angular Velocity
    
    % Parse output
    y(1) = d*p.omega;
    y(2) = 1-ga;
    y(3) = wf;
    y(4) = TbCum/n;
end




















