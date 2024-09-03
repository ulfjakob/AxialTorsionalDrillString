%%
clc
clear;
% close all;

p = parameters();


%% Sim stuff
T0 = 0;
T1 = 25;
p.t = T0:p.dt:T1;
p.Nt = numel(p.t);

% Init vals
% v0 = ones(p.P,1)*p.V0;
% v0(1) = 1;
% v0 = sin(linspace(0,2*pi,p.P).');
% v0(50:51) = [1; 1];
% h0 = linspace(p.Wf/(p.Ap*p.E) + p.L*p.ka/p.E*p.V0*p.rho, ...
%     p.Wf/(p.Ap*p.E),p.P).';
% h0 = linspace(2.8E-4,2E-4,p.P).';

v0 = ones(p.P,1)*p.V0 *0;
h0 = linspace(2.8E-4,2E-4,p.P).'*0;
p.K_a = 0;
p.Wf = 0;
p.V0 = 0;

% l0 = zeros(p.Pl,1);
l0 = linspace(0,1,p.Pl)*(2*pi)*p.V0/(p.omega0*p.N);
% Parse init states
x0(1:p.P)               = h0;        % Axial Strain
x0(p.P+1:2*p.P)         = v0;        % Axial Velocity
x0(2*p.P+1:2*p.P+p.Pl)  = l0;        % Depth of cut
x = zeros(2*p.P+p.Pl,p.Nt);
y = zeros(3,p.Nt);
x(:,1) = x0;

u.wb = 0;
%% Plot
[hAx,hBx] = setFigs(p);
%%
u.wb = 0;
u.tb = 0;
u.obEst = p.omega0;
for k=1:p.Nt

%     if p.t(k) > 1
        u.wb = 1;
%     end

    [x(:,k+1),y(1:3,k+1)] = waveStep_transportBRI(p,x(:,k),u);
    u.obEst = x(p.P,k+1);
    if mod(k,300)== 0
        set(hAx(1),'XData',x(1:p.P,k),'YData',p.x);
        set(hAx(2),'XData',x(p.P+1:2*p.P,k),'YData',p.x);

        set(hAx(3),'XData',linspace(0,1,p.Pl+1), ...
            'YData',[0;x(2*p.P+1:2*p.P+p.Pl,k)]);

        set(hBx(1),'XData',p.t,'YData',x(2*p.P,:)); % bit velocity
        
        set(hBx(2),'XData',p.t,'YData',x(2*p.P,:)); % bit velocity
        set(hBx(3),'XData',p.t,'YData',(y(1,:))/p.tN/p.N); % Norm depth of cut
        set(hBx(4),'XData',p.t,'YData',y(2,:)); % g nonlinearity
        drawnow
    %     pause(.1)
    %     F(k) = getframe(hFig);
    end
end

%%
% Axial limit cycle period
t_a = .2;

figure(10);
clf;
plot(p.t,x(2*p.P,1:p.Nt))
hold all
plot(p.t,(y(1,1:p.Nt)-p.dNom*p.N)/t_a)
plot(p.t,y(2,1:p.Nt))
l = legend('$v_b(t)$','$\frac{d(t)-\bar{d}}{t_a}$','$1-g_a(t)$');
set(l,'interpreter','latex');
title('Axial limit cycle')
xlabel('time (s)');
ylabel('');
xlim([83 85])

%%
zeta_a = p.Ap*p.E/p.c_a;
1/zeta_a
t_a = p.L/p.c_a
eta_a = exp(-2*p.ka*t_a)
eta_a/zeta_a

%%
figure(13);
clf;
plot(p.t,x(2*p.P,1:p.Nt))
hold on;
plot(p.t,1/zeta_a*exp(-p.ka*p.t),'--k');
% plot(p.t,p.t*0 +1/zeta_a,':k');
plot([t_a t_a t_a*2 t_a*2 t_a*3 t_a*3 t_a*4 t_a*4 t_a*5 t_a*5 t_a*6 t_a*6]*2,...
    [-1 1 1 -1 -1 1 1 -1 -1 1 1 -1],':k');
% plot(p.t,p.t*0 -1/zeta_a,':k');
plot(p.t,-1/zeta_a*exp(-p.ka*p.t),'--k');
plot(p.t,p.t*0,'-k');
l = legend('$v_b(t)$','$\pm \frac{e^{-k_at}}{\zeta_a}$',...
    '$2t_a,4t_a,6t_a,8t_a,\dots$');
set(l,'interpreter','latex');
title('Bit velocity step response')
xlabel('time (s)');
ylabel('$v_b$ (m/s)');
xlim([0 5])
ylim([-8e-6 8e-6]);











