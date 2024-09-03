%%
clc
clear;
close all;

%% Set parameter sets
tN_vec    = 1./[1 2 4 1 2 4 1 2 4];
eta_vec   = ones(size(tN_vec))*0.7;
Ka_vec = [10 10 10 20 20 20 40 40 40];
% MGVec = gainMargPipe(tN_vec,eta_vec,Ka_vec,'Z')
%%
% tN_vec    = .5;
% eta_vec   = -ones(size(tN_vec))*0.7;
% Ka_vec = 20;
%%
for i=1:numel(tN_vec)
    tN  = tN_vec(i)
    Ka  = Ka_vec(i)
    eta = eta_vec(i)
    [MG,wu] = ustabModes(tN,eta,Ka,'k');
    round(MG,3,'significant')
    wu = round(wu/(2*pi),3,'significant')
    keyboard
end

%%
% figure(10);
% TbLp = zeros(size(t));
% for k=2:numel(t)
%     TbLp(k) = TbLp(k-1)*.98 + y(4,k-1)*.02;
% end
% 
% clf;
% plot(t,y(5,:))
% hold on
% plot(t,y(1,:))
% plot(t,y(2,:)*p.Wf)
% plot(t,y(4,:))
% plot(t,TbLp(k),'-k','linewidth',2)
% plot(t,t*0+1,'--k')
% l = legend('$v_b(t)$','$d(t)$','$-\tilde{W}_f(t)$','$\tilde{T}_b(t)$',...
%     '$\tilde{T}_b(t) Lowpass filter$',...
%     'location','best');
% set(l,'interpreter','latex');
% title('Axial limit cycle')
% l = xlabel('Nondimensional time $\bar{t}$'); set(l,'interpreter','latex');
% l = ylabel('$v_b$ (m/s)'); set(l,'interpreter','latex');
% ylabel('');

%%
% h = figure(11);
% pdfmatlabfrag2(h,'LC33');
