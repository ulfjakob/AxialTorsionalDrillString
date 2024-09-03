function y = runOmegaSweep(p)


omegaVec = 10*[1./linspace(10,1,ceil(p.Nt/2)) 1./linspace(1,10,floor(p.Nt/2))];
p.omega = omegaVec(1);
p.tN = 1/p.omega;

% Init vals
% v0 = ones(p.P,1)*p.V0;
% v0(1) = 1;
% v0 = sin(linspace(0,2*pi,p.P).');
% v0(50:51) = [1; 1];
% h0 = linspace(p.Wf/(p.Ap*p.E) + p.L*p.ka/p.E*p.V0*p.rho, ...
%     p.Wf/(p.Ap*p.E),p.P).';
% h0 = linspace(2.8E-4,2E-4,p.P).';
v0 = ones(p.P,1)*p.v0;
h0 = linspace(p.k*p.v0,0,p.P).' +p.Wf+p.K_a*p.tN*p.v0;

% l0 = zeros(p.Pl,1);
l0 = linspace(0,1,p.Pl)/(p.omega);
% Parse init states
x0(1:p.P)               = h0;        % Axial Strain
x0(p.P+1:2*p.P)         = v0;        % Axial Velocity
x0(2*p.P+1:2*p.P+p.Pl)  = l0;        % Depth of cut
x = zeros(2*p.P+p.Pl,p.Nt);
y = zeros(4,p.Nt);
x(:,1) = x0;

u.wb = 0;
%% Plot
[hAx,hBx] = setFigs(p);
%%
TbLp = nan(1,p.Nt);
TbLp(1) = 0;
u.wb = 0;
u.tb = 0;
for k=1:p.Nt
    
    if p.t(k)>=1
        p.v0 = 1;
    end
    
    p.omega = omegaVec(k);
    p.tN = 1/p.omega;
    
    [x(:,k+1),y(1:4,k+1)] = waveStep_transportBRI(p,x(:,k),u);
    u.obEst = x(p.P,k+1);
    
    if k>1; TbLp(k) = TbLp(k-1)*.95 + y(4,k)*.05; end
    
    if mod(k,p.drawPeriod)== 0
        set(hAx(1),'XData',x(1:p.P,k),'YData',p.x);
        set(hAx(2),'XData',x(p.P+1:2*p.P,k),'YData',p.x);

        set(hAx(3),'XData',linspace(0,1,p.Pl+1), ...
            'YData',[0;x(2*p.P+1:2*p.P+p.Pl,k)]);

        set(hBx(1),'XData',p.t,'YData',x(2*p.P,:)); % bit velocity
        set(hBx(2),'XData',p.t,'YData',(y(1,:))); % Norm depth of cut
        set(hBx(3),'XData',p.t,'YData',(y(2,:))*p.Wf); % wf force
        set(hBx(4),'XData',p.t,'YData',(y(4,:))); % Tb
        set(hBx(5),'XData',p.t,'YData',TbLp); % Tb Lp

        drawnow
    %     pause(.1)
    %     F(k) = getframe(hFig);
    end
end

%% Log vb
y(5,1:p.Nt) = x(2*p.P,1:p.Nt);
y(6,1:p.Nt) = omegaVec;
%%
% Axial limit cycle period
% figure(10);
figure;
clf;
plot(p.t,x(2*p.P,1:p.Nt))
hold all
plot(p.t,(y(1,1:p.Nt)))
plot(p.t,y(2,1:p.Nt)*p.Wf)
plot(p.t,y(4,1:p.Nt))
plot(p.t,p.t*0+1,'--k')
l = legend('$v_b(t)$','$d(t)$','$-\tilde{W}_f(t)$','$\tilde{T}_b(t)$',...
    'location','best');
set(l,'interpreter','latex');
title('Axial limit cycle')
l = xlabel('Nondimensional time $\bar{t}$'); set(l,'interpreter','latex');
l = ylabel('$v_b$ (m/s)'); set(l,'interpreter','latex');
ylabel('');
xlim([55 60])

%%
% gi = y(2,:);
% k = 1;
% stick = 0;
% i1 = 0;
% TbCum = 0;
% for i=1:numel(gi)
%     if stick == 0 && gi(i)
%         TbAvg(k) = TbCum/((i-i1)*p.dt);
%         k = k+1;
%         i1 = i;
%         stick = 1;
%         TbCum = y(4,i);
%     elseif stick == 1 && gi(i) > 0
%         TbCum = TbCum + y(4,i);
%     elseif stick == 1 && gi(i) == 0
%         stick = 0;
%     end
% end
% %
% figure(41)
% plot(TbAvg)
% %
% try
%     TbAvg = mean(TbAvg(end-10:end));
% catch
%     TbAvg = nan;
% end
%%
% figure(13);
% clf;
% plot(p.t,-1+x(2*p.P,1:p.Nt))
% hold on;
% plot(p.t,exp(-p.k*p.t),'--k');
% plot([1 1 2 2 3 3 4 4 5 5 6 6]*2 -1,...
%     [-1 1 1 -1 -1 1 1 -1 -1 1 1 -1],':k');
% plot(p.t,-exp(-p.k*p.t),'--k');
% plot(p.t,p.t*0,'-k');
% l = legend('$v_b(t)$','$\pm \frac{e^{-k_at}}{\zeta_a}$',...
%     '$1t_a,3t_a,5t_a,7t_a,\dots$');
% set(l,'interpreter','latex');
% title('Bit velocity step response')
% l = xlabel('Nondimensional time $\bar{t}$'); set(l,'interpreter','latex');
% l = ylabel('$v_b$ (m/s)'); set(l,'interpreter','latex');
% xlim([0 20])

% h = figure(13);
% pdfmatlabfrag2(h,'stepResp');

%%
% figure(11);
% clf;
% plot(p.t,x(2*p.P,1:p.Nt))
% hold all
% plot(p.t,(y(1,1:p.Nt)-p.dNom*p.N)/t_a)
% plot(p.t,y(2,1:p.Nt))
% l = legend('$v_b(t)$','$\frac{d(t)-\bar{d}}{t_a}$','$1-g_a(t)$');
% set(l,'interpreter','latex');
% title('Axial limit cycle')
% xlabel('time (s)');
% ylabel('');
% xlim([20 21])
% %%
% figure(12);
% clf;
% plot(p.t,x(2*p.P,1:p.Nt))
% hold all
% plot(p.t,(y(1,1:p.Nt)-p.dNom*p.N)/t_a)
% plot(p.t,y(2,1:p.Nt))
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
















