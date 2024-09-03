function [hAx,hBx,hFig] = setFigs(p)

hFig(1) = figure(1);
clf;
subplot(131);
hAx(1) = plot(nan(size(p.x)),p.x);
% xlim([0.90 2.00]*1e-4);
ylim([0 p.L]);
title('Axial Strain'); xlabel('(m/m)')
ylabel('MD (m)');
set (gca,'Ydir','reverse')

subplot(132); 
hAx(2) = plot(nan(size(p.x)),p.x);
ylim([0 p.L]);
xlim([p.V0-.05 5*p.V0+.05]);
title('Axial Velocity'); xlabel('(m/s)')
ylabel('MD (m)');
set (gca,'Ydir','reverse')

subplot(133);
hAx(3) = plot(linspace(0,1,p.Pl),nan(p.Pl,1));
hold on;
dNom = (2*pi)*p.V0/(p.omega0*p.N);
plot([.9 1],dNom*[1 1]);
xlim([0 1]); 
% ylim([-dNom*.2 dNom*2]);
title('\lambda'); xlabel('\theta')
ylabel('Depth of cut (m)');

%%%%%%%%%%%%%%%%%%%%
hFig(2) = figure(2);
clf;
subplot(211);
hBx(1) = plot(nan(size(p.x)),p.x);
title('bit axial velocity');

subplot(212);
plot(p.t, ones(size(p.t))*p.V0,'--k' );
hold on;
hBx(2) = plot(nan(size(p.t)),p.t);
hBx(3) = plot(nan(size(p.t)),p.t);
hBx(4) = plot(nan(size(p.t)),p.t);
% ylim([-.2 2]*p.V0)
title('Depth of cut');
legend('Nominal','Bit velocity','Depth of cut','g nonlin');





