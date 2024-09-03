function [hAx,hBx,hFig] = setFigs(p)

hFig = figure(1);
clf;
subplot(231);
hAx(1) = plot(nan(size(p.x)),p.x);
xlim([0.90 2.00]*1e-4);
ylim([0 p.L]);
title('Axial Strain'); xlabel('(m/m)')
ylabel('MD (m)');
set (gca,'Ydir','reverse')

subplot(232); 
hAx(2) = plot(nan(size(p.x)),p.x);
ylim([0 p.L]);
xlim([p.V0-.05 5*p.V0+.05]);
title('Axial Velocity'); xlabel('(m/s)')
ylabel('MD (m)');
set (gca,'Ydir','reverse')

subplot(234);
hAx(3) = plot(nan(size(p.x)),p.x);
ylim([0 p.L]); 
xlim([-.01 .03]); 
title('Angular Strain'); xlabel('(rad/m)')
ylabel('MD (m)');
set (gca,'Ydir','reverse')

subplot(235); 
hAx(4) = plot(nan(size(p.x)),p.x);
ylim([0 p.L]); 
xlim([-1*p.omega0-5 p.omega0*4+5]);
title('Angular Velocity'); xlabel('(rad/s)')
ylabel('MD (m)');
set (gca,'Ydir','reverse')

subplot(233);
hAx(5) = plot(linspace(0,1,p.Pl),nan(p.Pl,1));
hold on;
dNom = (2*pi)*p.V0/(p.omega0*p.N);
plot([.9 1],dNom*[1 1]);
xlim([0 1]); 
ylim([-dNom*.2 dNom*2]);
title('\lambda'); xlabel('\theta')
ylabel('Depth of cut (m)');

%%%%%%%%%%%%%%%%%%%%
figure(2);
clf;
subplot(221);
hBx(1) = plot(nan(size(p.x)),p.x);
title('bit axial velocity');

subplot(222); 
hBx(2) = plot(nan(size(p.x)),p.x);
title('Bit angular velocity');

subplot(212);
hBx(3) = plot(nan(size(p.x)),p.x);
hold on;
hBx(4) = plot(nan(size(p.x)),p.x);
hBx(5) = plot(nan(size(p.x)),p.x);
title('Depth of cut');





