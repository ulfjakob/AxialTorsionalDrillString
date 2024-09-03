function VbAvg_contour(omMat,W0Mat,V0AvgMat,p)


% contourf(omMat,W0Mat,V0AvgMat,[0 .1:.2:10 .99])
% contourf(omMat,W0Mat,V0AvgMat,'ShowText','on')

% contourVec = [0 logspace(-1,3,13)];
contourVec = [0 0.3 0.5 0.7 0.9 1 1.05 1.1 1.15  1.2 1.25 1.3];
contourf(omMat,W0Mat,V0AvgMat,'LevelList',contourVec,'ShowText','on')
caxis([0.25,1.20])
set(gca,'yscale','log')
set(gca,'xscale','log')

% ylim([.1 10^1.5])

%%
ylabel('W_0-W_f');
xlabel('\Omega_0=1/t_N');

% tString = sprintf('$K_a=%d, \\ \\eta_a=%.2f \\ \\eta_t=%.2f$',p.K_a,p.eta_a,...
%     p.eta_t);
% l = title(tString); set(l,'interpreter','latex');

% setLatexLabels2

