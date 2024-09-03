function [hAx,hBx,hFig] = setFigs(p)
v0 = max(p.v0,1);

hFig = figure(1);
clf;
subplot(231);
hAx(1) = plot(nan(size(p.x)),p.x);
d0  = p.tN0*v0;
NomStrain = p.Wf+p.K_a*d0+[0-1,p.ka*v0*2+1];
xlim(NomStrain);
ylim([0 1]);
title('Axial Strain'); xlabel('(m/m)')
ylabel('MD');
set (gca,'Ydir','reverse')

subplot(232); 
hAx(2) = plot(nan(size(p.x)),p.x);
xlim([-1 2]);
ylim([0 1]);
title('Axial Velocity'); xlabel('(m/s)')
ylabel('MD (m)');
set (gca,'Ydir','reverse')

subplot(234);
hAx(3) = plot(nan(size(p.x)),p.x);
omega0 = 1/p.tN0;
% NomAngStrain = p.Wf*p.beta+p.K_t*d0+[0-1,p.kt*p.tt^2*omega0*2+1];
NomAngStrain = p.Wf*p.beta+d0+[-d0,p.kt*omega0*2+d0];
xlim(NomAngStrain);
ylim([0 1]); 
title('Angular Strain'); xlabel('(rad/m)')
ylabel('MD (m)');
set (gca,'Ydir','reverse')

subplot(235); 
hAx(4) = plot(nan(size(p.x)),p.x);
ylim([0 1]); 
xlim([-1 2]*omega0);
title('Angular Velocity'); xlabel('(rad/s)')
ylabel('MD (m)');
set (gca,'Ydir','reverse')

subplot(233);
hAx(5) = plot(linspace(0,1,p.Pl),nan(p.Pl,1));
hold on;
dNom = v0*p.tN0;
plot([.9 1],dNom*[1 1]);
plot([0 1],[0 0],'--k');
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





