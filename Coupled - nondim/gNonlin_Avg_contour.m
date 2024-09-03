function gNonlin_Avg_contour(omMat,W0Mat,gBar,p)


% contourf(omMat,W0Mat,V0AvgMat,[0 .1:.2:10 .99])
% contourf(omMat,W0Mat,V0AvgMat,'ShowText','on')

% contourVec = [-2 -1 -.5 -.1 -.01 0 logspace(-1,3,13)];
% contourVec = [.5 1 1.1 1.2 1.5 2 4 6];
contourVec = [-10 -7 -4 -1 0 1 4 7 10];
% contourVec = [0 .8 .9 .98 1 1.1 linspace(1.5,3,11)];
contourf(omMat,W0Mat,gBar,contourVec,'ShowText','on')
set(gca,'yscale','log')
set(gca,'xscale','log')

% ylim([.1 10^1.5])

%%
l = ylabel('$W_0-W_f$'); set(l,'interpreter','latex');
l = xlabel('$\Omega_0=1/t_N$'); set(l,'interpreter','latex');

% tString = sprintf('$K_a=%d, \\ \\eta_a=%.2f \\ \\eta_t=%.2f$',p.K_a,p.eta_a,...
%     p.eta_t);
% l = title(tString); set(l,'interpreter','latex');

% setLatexLabels2

