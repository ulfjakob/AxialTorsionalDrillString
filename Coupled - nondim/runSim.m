function [y,t] = runSim(p,plottings)

if nargin == 1
    plottings = 1;
end

d0  = p.tN0*p.v0;
% Init vals
xPos = linspace(0,1,p.P).';

v0 = xPos*0 + p.v0;
omega0 = 1/p.tN0;
if p.eta_a<0
    h0 = p.w0 - flipud(linspace(p.ka*p.v0,0,p.P).');
else
    h0 = linspace(p.ka*p.v0,0,p.P).' +p.Wf+p.K_a/omega0*p.v0;
end
l0 = linspace(0,1,p.Pl)*p.v0/omega0;
o0 = xPos*0 + omega0;
f0 = d0 + p.Wf*p.beta + p.kt*omega0*(1-xPos);

% Parse init states
x0(1:p.P)               = h0;        % Axial Strain
x0(p.P+1:2*p.P)         = v0;        % Axial Velocity
x0(2*p.P+1:3*p.P)       = f0;        % Angular Strain
x0(3*p.P+1:4*p.P)       = o0;        % Angular Velocity
x0(4*p.P+1:4*p.P+p.Pl)  = l0;        % Depth of cut
x = zeros(4*p.P+p.Pl,p.Nt);
y = zeros(4,p.Nt);
x(:,1) = x0;

%% Plot
if plottings
    [hAx,hBx] = setFigs(p);
end
%%
u.wb = 0;
u.tb = 0;
for k=1:p.Nt-1
    
    if p.t(k)<=15
        if p.eta_a<0
%             p.w0 = p.Wf +(p.k + p.K_a*1 )*p.v0;
            u.wb = 1;
        else
%             p.v0 = 1;
            u.wb = 1;
        end
    else
        u.wb = 0;
    end
    
%     [x(:,k+1),y(1:3,k+1)] = waveStepLFSG(p,x(:,k),u);
    [x(:,k+1),y(1:3,k+1)] = waveStepUpwind(p,x(:,k),u);
    u.obEst = x(p.P,k+1);
    
    if plottings && mod(k,p.drawPeriod) == 0
        set(hAx(1),'XData',x(1:p.P,k)           ,'YData',p.x);
        set(hAx(2),'XData',x(p.P+1:2*p.P,k)     ,'YData',p.x);
        set(hAx(3),'XData',x(2*p.P+1:3*p.P,k)   ,'YData',p.x);
        set(hAx(4),'XData',x(3*p.P+1:4*p.P,k)   ,'YData',p.x);

        set(hAx(5),'XData',linspace(0,1,p.Pl+1), ...
            'YData',[0;x(4*p.P+1:4*p.P+p.Pl,k)]);

        %
        set(hBx(1),'XData',p.t,'YData',x(2*p.P,:));     % bit velocity
        set(hBx(2),'XData',p.t,'YData',x(4*p.P,:));     % bit angular velocity
        set(hBx(3),'XData',p.t,'YData',x(2*p.P,:));
        set(hBx(4),'XData',p.t,'YData',(y(1,:)));       % Norm depth of cut
        set(hBx(5),'XData',p.t,'YData',y(2,:)*p.Wf);    % wf force

        drawnow
    %     F(k) = getframe(hFig);
    end
end

%% Log vb
y(4,1:p.Nt) = x(2*p.P,1:p.Nt);
y(5,1:p.Nt) = x(4*p.P,1:p.Nt);
t = p.t(1:p.Nt);

