function surfMargins(X,Y,Z)

% tString = sprintf('$\\bar{\\zeta}$ = %1.2f, $\\eta$=%1.2f',eta);
tString = '';
figure; clf;
% contourf(tpVec,t_nVec,db(Z),[-6:.5:3]*10,'ShowText','on')
% colormap default
surf(X,Y,db(abs(Z)+1e-2))
% caxis([.1,1.5])
caxis([-30,10])
l = title(tString);
set(l,'interpreter','latex');
set(gca,'yscale','log')
l = xlabel('$\eta_a$'); set(l,'interpreter','latex');
l = ylabel('$\bar{t}_N$'); set(l,'interpreter','latex');
l = zlabel('$M_{G,a}(dB)$'); set(l,'interpreter','latex');
setLatexLabelsSurf('y')

ylim([min(Y) max(Y)]);

end


function setLatexLabelsSurf(axis)

if nargin == 0
    axis = 'b';
end
if strcmp(axis,'y') || strcmp(axis,'b')
    x=xlim;
    xp=x(2)+abs(diff(xlim))*0.03;   % For other side on 3d plots
    y=get(gca,'ytick');
    z = zlim;
    set(gca,'yticklabel',['$$'],'xtickmode','manual');
    for j=1:length(y)
      h=text(xp,y(j),z(1),['$10^{' num2str(log10(y(j))) '}$']);
      set(h,'HorizontalAlignment','center','fontsize',get(gca,'fontsize'))
    end
end

if strcmp(axis,'x') || strcmp(axis,'b')
    y=ylim;
    yp=y(1)-abs( log(diff(ylim))*0.008 );
    % yp = .08;
    x=get(gca,'xtick');
    set(gca,'xticklabel',['$$'],'xtickmode','manual');
    for j=1:length(x)
      h=text(x(j),yp,['$10^{' num2str(log10(x(j))) '}$']);
      set(h,'HorizontalAlignment','center','fontsize',get(gca,'fontsize'))
    end
end
end