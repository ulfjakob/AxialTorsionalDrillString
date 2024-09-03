function [y,t] = runSim(p)

% Init vals
v0 = ones(p.P,1)*p.v0;
if p.eta<0
    h0 = p.w0 - flipud(linspace(p.k*p.v0,0,p.P).');
    l0 = linspace(0,1,p.Pl)*p.v0;
else
    h0 = linspace(p.k*p.v0,0,p.P).' +p.Wf+p.K_a/p.omega*p.v0;
    l0 = linspace(0,1,p.Pl)*p.v0/(p.omega);
end

% l0 = zeros(p.Pl,1);

% Parse init states
x0(1:p.P)               = h0;        % Axial Strain
x0(p.P+1:2*p.P)         = v0;        % Axial Velocity
x0(2*p.P+1:2*p.P+p.Pl)  = l0;        % Depth of cut
x = zeros(2*p.P+p.Pl,p.Nt);
y = zeros(4,p.Nt);
x(:,1) = x0;

%% Plot
if ~(p.drawPeriod==inf)
    [hAx,hBx] = setFigs(p);
end
%%
TbLp = nan(1,p.Nt);
TbLp(1) = 0;
% u.wb = 0;
% u.wb = 1;
% u.tb = 0;
for k=1:p.Nt-1
    
    if p.t(k)<=5
        if p.eta<0
%             p.w0 = p.Wf +(p.k + p.K_a*1 )*p.v0;
            u.wb = 2;
        else
            p.v0 = 1;
            u.wb = 2;
        end
    else
        u.wb = 0;
    end
%     if p.t(k)>3
%         p.v0 = 1;
%     end
    
    [x(:,k+1),y(1:4,k+1)] = waveStep_transportBRI(p,x(:,k),u);
    u.obEst = x(p.P,k+1);
    
    if k>1; TbLp(k) = TbLp(k-1)*.95 + y(4,k)*.05; end
    
    if mod(k,p.drawPeriod)== 0
        set(hAx(1),'XData',x(1:p.P,k),'YData',p.x);
        set(hAx(2),'XData',x(p.P+1:2*p.P,k),'YData',p.x);

        set(hAx(3),'XData',linspace(0,1,p.Pl+1), ...
            'YData',[0;x(2*p.P+1:2*p.P+p.Pl,k)]);

        set(hBx(1),'XData',p.t,'YData',x(2*p.P,:)); % bit velocity
        set(hBx(2),'XData',p.t,'YData',(y(1,:))); % Norm depth of cut
        set(hBx(3),'XData',p.t,'YData',(y(2,:))*p.Wf); % wf force
%         set(hBx(4),'XData',p.t,'YData',(y(4,:))); % Tb
%         set(hBx(5),'XData',p.t,'YData',TbLp); % Tb Lp

        drawnow
    %     pause(.1)
    %     F(k) = getframe(hFig);
    end
end

%% Log vb
y(5,1:p.Nt) = x(2*p.P,1:p.Nt);
t = p.t(1:p.Nt);


