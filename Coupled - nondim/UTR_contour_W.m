function UTR_contour_W(omMat,W0Mat,UTR_Mat,p)

UTR_Mat_c = min(UTR_Mat,1);
UTR_Mat_c = UTR_Mat_c+(UTR_Mat_c>=1)*.5;
% UTR_Mat_c = UTR_Mat;
contourf(omMat,W0Mat,UTR_Mat_c,[0 .1:.2:.9 .99])
% contourf(omMat,KtMat,UTR_Mat_c,'ShowText','on')
set(gca,'yscale','log')
set(gca,'xscale','log')

l = ylabel('W_0-W_f');
% set(l,'interpreter','latex');


l = xlabel('\Omega_0=1/t_N');
% set(l,'interpreter','latex');

% tString = sprintf('$K_a=%d, \\ \\eta_a=%.2f \\ \\eta_t=%.2f$',p.K_a,p.eta_a,...
%     p.eta_t);
% l = title(tString); set(l,'interpreter','latex');
%%
x1t1 = 0.3;
y1t1 = 150.0;
x1t2 = 0.08;
y1t2 = 4;
rot = 15;


txt1 = {'Stick slip limit cycle'};
l = text(x1t1,y1t1,txt1); 
set(l,'fontSize',16,'color',[.8 .0 .0],'rotation',rot);
txt2 = {'Stable torsional dynamics'};
l = text(x1t2,y1t2,txt2);
set(l,'fontSize',16,'color',[.8 .0 .0],'rotation',rot);
%%
% setLatexLabels2

