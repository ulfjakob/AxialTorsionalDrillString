function UTR_contour(omMat,V0Mat,UTR_Mat,p)

UTR_Mat_c = min(UTR_Mat,1);
UTR_Mat_c = UTR_Mat_c+(UTR_Mat_c>=1)*.5;
% UTR_Mat_c = UTR_Mat;
contourf(omMat,V0Mat,UTR_Mat_c,[0 .1:.2:.9 .99])
% contourf(omMat,KtMat,UTR_Mat_c,'ShowText','on')
set(gca,'yscale','log')
set(gca,'xscale','log')

if p.eta_a>0
    l = ylabel('$V_0$'); set(l,'interpreter','latex');
else
    l = ylabel('Average $V_b$'); set(l,'interpreter','latex');
end

l = xlabel('$\Omega_0=1/t_N$'); set(l,'interpreter','latex');

% tString = sprintf('$K_a=%d, \\ \\eta_a=%.2f \\ \\eta_t=%.2f$',p.K_a,p.eta_a,...
%     p.eta_t);
% l = title(tString); set(l,'interpreter','latex');

x1t1 = .26;
y1t1 = 1.0;
x1t2 = 2.23;
y1t2 = 0.2;
rot = 55;

txt1 = {'Stick slip limit cycle'};
l = text(x1t1,y1t1,txt1); 
set(l,'interpreter','latex','fontSize',16,'color',[.8 .0 .0],'rotation',rot);
txt2 = {'Stable torsional dynamics'};
l = text(x1t2,y1t2,txt2);
set(l,'interpreter','latex','fontSize',16,'color',[.8 .0 .0],'rotation',rot);

setLatexLabels2

