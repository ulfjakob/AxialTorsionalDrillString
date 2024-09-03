%%
clc
clear;
% close all;

p = parameters();
clear waveStep
%% Sim stuff
T0 = 0;
T1 = 25;
t = T0:p.dt:T1;
Nt = numel(t);

% Init vals
v0 = ones(p.P,1)*p.V0;
% v0(1) = 1;
% v0 = sin(linspace(0,2*pi,p.P).');
% v0(50:51) = [1; 1];
h0 = linspace(p.Wf/(p.Ap*p.E)*1.8 + p.L*p.ka/p.E*p.V0*p.rho, ...
    p.Wf/(p.Ap*p.E),p.P).';
o0 = ones(p.P,1)*p.omega0;
f0 = linspace(p.Tf*p.Wf/(p.Jp*p.G) + p.L*p.kt/p.G*p.omega0*p.rho, ...
    p.Tf*p.Wf/(p.Jp*p.G),p.P).';
% l0 = zeros(p.Pl,1);
l0 = linspace(0,1,p.Pl)*(2*pi)*p.V0/(p.omega0*p.N);
% Parse init states
x0(1:p.P)               = h0;        % Axial Strain
x0(p.P+1:2*p.P)         = v0;        % Axial Velocity
x0(2*p.P+1:3*p.P)       = f0;        % Angular Strain
x0(3*p.P+1:4*p.P)       = o0;        % Angular Velocity
x0(4*p.P+1:4*p.P+p.Pl)  = l0;        % Depth of cut
x = zeros(4*p.P+p.Pl,Nt);
y = zeros(3,Nt);
x(:,1) = x0;

%% Plot
[hAx,hBx] = setFigs(p);
%%
u.wb = 0;
u.tb = 0;
u.obEst = p.omega0;
for k=1:Nt

%     if ismember(k,1000+[0])
% %     if 0
%         W_impulse = 2e1 ;
%         u.wb = W_impulse*10;
%         u.tb = W_impulse*10;
%     else
% %         u.wb = 0;
% %         u.tb = 0;
%     end

    [x(:,k+1),y(1:3,k+1)] = waveStep_transportBRI(p,x(:,k),u);
    u.obEst = x(p.P,k+1);
    if mod(k,1)== 0
        set(hAx(1),'XData',x(1:p.P,k),'YData',p.x);
        set(hAx(2),'XData',x(p.P+1:2*p.P,k),'YData',p.x);
        set(hAx(3),'XData',x(2*p.P+1:3*p.P,k),'YData',p.x);
        set(hAx(4),'XData',x(3*p.P+1:4*p.P,k),'YData',p.x);

        set(hAx(5),'XData',linspace(0,1,p.Pl+1), ...
            'YData',[0;x(4*p.P+1:4*p.P+p.Pl,k)]);

        set(hBx(1),'XData',t,'YData',x(2*p.P,:)); % bit velocity
        set(hBx(2),'XData',t,'YData',x(4*p.P,:)); % bit angular velocity
        set(hBx(3),'XData',t,'YData',x(2*p.P,:));
        set(hBx(4),'XData',t,'YData',(y(1,:)-p.dNom*p.N));
        set(hBx(5),'XData',t,'YData',y(2,:));
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
plot(t,x(2*p.P,1:Nt))
hold all
plot(t,(y(1,1:Nt)-p.dNom*p.N)/t_a)
plot(t,y(2,1:Nt))
l = legend('$v_b(t)$','$\frac{d(t)-\bar{d}}{t_a}$','$1-g_a(t)$');
set(l,'interpreter','latex');
title('Axial limit cycle')
xlabel('time (s)');
ylabel('');
xlim([24 25])

% h = figure(10);
% pdfmatlabfrag2(h,'axialLS_30RPM');

%%
% figure(11);
% clf;
% plot(t,x(2*p.P,1:Nt))
% hold all
% plot(t,(y(1,1:Nt)-p.dNom*p.N)/t_a)
% plot(t,y(2,1:Nt))
% l = legend('$v_b(t)$','$\frac{d(t)-\bar{d}}{t_a}$','$1-g_a(t)$');
% set(l,'interpreter','latex');
% title('Axial limit cycle')
% xlabel('time (s)');
% ylabel('');
% xlim([20 21])
% %%
% figure(12);
% clf;
% plot(t,x(2*p.P,1:Nt))
% hold all
% plot(t,(y(1,1:Nt)-p.dNom*p.N)/t_a)
% plot(t,y(2,1:Nt))
% l = legend('$v_b(t)$','$\frac{d(t)-\bar{d}}{t_a}$','$1-g_a(t)$');
% set(l,'interpreter','latex');
% title('Axial limit cycle')
% xlabel('time (s)');
% ylabel('');
% xlim([18 19])
%% Save video
% v = VideoWriter('ss.avi');
% open(v);
% writeVideo(v,F);
% close(v)

%%
% h = figure(2);
% pdfmatlabfrag2(h,'hh');
% h = figure(3);
% pdfmatlabfrag2(h,'Ha');
% h = figure(4);
















