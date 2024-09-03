function [hAx,hBx,hFig] = setFigs(p)

hFig = figure(1);
clf;
subplot(131);
hAx(1) = plot(nan(size(p.x)),p.x);
ylim([0 1]);
v0 = 1;
xlim([-10 p.k*v0+10]+p.Wf+p.K_a/p.omega*v0);
title('Axial Strain'); xlabel('(m/m)')
ylabel('MD (m)');
set (gca,'Ydir','reverse')

subplot(132); 
hAx(2) = plot(nan(size(p.x)),p.x);
ylim([0 1]);
xlim([0 2]);
title('Axial Velocity'); xlabel('(m/s)')
ylabel('MD (m)');
set (gca,'Ydir','reverse')

subplot(133);
hAx(3) = plot(linspace(0,1,p.Pl),nan(p.Pl,1));
hold all;
dNom = 1/p.omega;
plot([.9 1],dNom*[1 1]);
xlim([0 1]); 
% ylim([-dNom*.2 dNom*2]);
title('\lambda'); xlabel('\theta')
ylabel('Depth of cut (m)');

%%%%%%%%%%%%%%%%%%%%
figure(2);
clf;
plot(p.t, ones(size(p.t)),'--k' );
hold all;
hBx(1) = plot(nan(size(p.t)),p.t);
hBx(2) = plot(nan(size(p.t)),p.t);
hBx(3) = plot(nan(size(p.t)),p.t);
hBx(4) = plot(nan(size(p.t)),p.t);
hBx(5) = plot(nan(size(p.t)),p.t,'--');
% ylim([-1 3]*1)
ylim([-10 20]*1)
title('Depth of cut');
legend('Nominal','Bit velocity','Depth of cut','$w_f$','$T_b$',...
    '$T_b Avg.$');




