%%
clc
clear;
close all;


%% Stability Map Location
% stabilityMap_pipe([],[]);

% figure(21);
% plot([1 1]*.7,[1 1/10],'-r','linewidth',2);
% 
% str1 = ['$\Omega_b$ sweep'];
% l = text(.65,.55,str1,'color','red');
% set(l,'interpreter','latex','rotation', 90);
% 
% h = figure(21);
% pdfmatlabfrag2(h,'stabMap');
%% Negative Gain Margins
% Set parameter sets
omega = 10*[1./linspace(10,1,300)];
tN_vec    = 1./omega;
eta_vec   = ones(size(tN_vec))*0.7;
Ka_vec    = ones(numel(tN_vec),1)*20;
MGVec = gainMargPipe(tN_vec,eta_vec,Ka_vec,'k');

figure(12);
semilogx(omega,MGVec);
%% Simulation analysis
p = parameters;

%% Sim stuff
T0 = 0;
T1 = 60*40;
p.t = T0:p.dt:T1;
p.Nt = numel(p.t);

i = 1;
p.tN = tN_vec(i);
p.omega = 1/(p.tN);
p.eta = eta_vec(i);
p.k = -log(p.eta)/2;
p.K_a = Ka_vec(i);

y = runOmegaSweep(p);
%%
% fp = sprintf('sweep%d',round(p.K_a));
% save(fp,'y')
%%
TbLp = zeros(size(p.t));
tLp = 0.010;
for k=2:numel(p.t)
    TbLp(k) = TbLp(k-1)*(1-tLp) + y(4,k-1)*tLp;
end
Tbip = zeros(size(p.t));
for i=1:numel(p.t)-1
    k = p.Nt-i;
    Tbip(k) = Tbip(k+1)*(1-tLp) + y(4,k+1)*tLp;
end
%%
figure(11);
clf
% subplot(211)
semilogx(y(6,1:ceil(p.Nt/2)),TbLp(1:ceil(p.Nt/2)));
hold on
semilogx(y(6,ceil(p.Nt/2)+1:end),TbLp(ceil(p.Nt/2):end));
% plot(y(6,ceil(p.Nt/2)+1:end),TbLp(ceil(p.Nt/2):end));
% plot(y(6,1:ceil(p.Nt/2)),Tbip(1:ceil(p.Nt/2)));
% plot(y(6,ceil(p.Nt/2)+1:end),Tbip(ceil(p.Nt/2):end));
semilogx(y(6,ceil(p.Nt/2)+1:end),1./y(6,ceil(p.Nt/2)+1:end),...
    'linewidth',2);
xlim([1 10]);
ylim([-1.5 1])
l = legend('Increasing $\Omega_b$','Decreasing $\Omega_b$',...
    '$g(V_b) = 1$','location','best');
set(l,'interpreter','latex');

setLatexLabels
xlabel('$\Omega_b$');
ylabel('Smoothed $T_b-\beta$');

% subplot(212)
% semilogx(omega,MGVec);
% ylim([5 30])
h = figure(11);
pdfmatlabfrag2(h,'sweep2');

%%
figure(10);
clf;
plot(p.t,y(5,2:end))
hold on;
plot(p.t,y(4,2:end))
plot(p.t,TbLp,'linewidth',2);
plot(p.t,1./y(6,1:p.Nt),'linewidth',1);
ylim([-5 5])
xlim([0 p.t(end)]);
l = legend('$v_b$','$T_b-\beta$','$T_b-\beta$ Smoothed',...
    '$t_N=\frac{V_b}{\Omega_b}$','location','best');
set(l,'interpreter','latex');
l = xlabel('Non-dimensional time'); set(l,'interpreter','latex')
l = ylabel('$T_b-\beta, \ V_b, \ t_N$'); set(l,'interpreter','latex')


% h = figure(10);
% pdfmatlabfrag2(h,'sweep1');

%%














