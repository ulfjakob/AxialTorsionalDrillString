function stabilityMap_pipe(tN_c,eta_c)

% nw = 1e3;
nw = 1e4;
w = linspace(3e-3,10*2*pi,nw)';
% Angular frequency the Laplace variable is evaluated at
% w = linspace(3e-3,15*2*pi,nw)';
s = 1j*w;

%% Physical Params
Ka = 1;
% etaVec = linspace(.03,.97,20);
% etaVec = linspace(-.97,.97,60);
etaVec = linspace(.03,.97,120);
t_nVec = logspace(-1,1,120);
% kaVec = logspace(-1,.5,5);
% t_nVec = logspace(-1,1,5);
%%
Z = nan(numel(t_nVec),numel(etaVec));
Z_wu = nan(numel(t_nVec),numel(etaVec));
for i = 1:numel(t_nVec)
    i
    t_n = t_nVec(i);
    for j = 1:numel(etaVec)
        eta = etaVec(j);
        Z_Lbar = (1+eta)/(1-eta);
        
        Zc  =       1;
        Ga   =      s;
        delay = (1-exp(-t_n*s));
        ga_i = (1+Z_Lbar*tanh(Ga))./(Z_Lbar+1*tanh(Ga));
        Ga_i = ga_i .*Ka./s.*delay;

        [unstab,dG,wu] = nyqStab(Ga_i,w);
        
        %%
%         db(dG)
%         wu/pi*2
%         figure(3); clf;
%         bodelin(Ga_i,w)
%         subplot(211)
%         ylim([-30 30])
        %%
        if isempty(dG)
            Z(i,j) = 0;
            Z_wu(i,j) = 0;
        else
            [Z(i,j)] = max(dG);
            if db(Z(i,j)) < -30
                Z_wu(i,j) = 0;
            else
                [Z_wu(i,j)] = wu;
            end
        end
    end
end
%%
% surfMargins(etaVec,t_nVec,Z)
%%
figure(21); clf;
contourf(etaVec,t_nVec,db(Z),[-4:1:3]*10,'ShowText','on')
hold all;
plot(eta_c,tN_c,'xr','linewidth',2)
plot(eta_c,tN_c,'or','linewidth',2)
l = title('Inverse gain margin (dB) $M_G = \max_{\varpi\in \varpi_{180}} \Big|G_a(j\varpi) \Big|$');
set(l,'interpreter','latex');
set(gca,'yscale','log')
l = xlabel('$\eta_a$'); set(l,'interpreter','latex');
l = ylabel('$t_N=1/\Omega_0$'); set(l,'interpreter','latex');
setLatexLabels2('y')
%%
figure(22); clf;
% contourf(etaVec,t_aVec,db(Z_wu),[-30 -20 -10 0 10 20],'ShowText','on')
contourf(etaVec,t_nVec,(Z_wu/pi*2)/2,[.1 0.8:1:12]);
hold all;
plot(eta_c,tN_c,'xr','linewidth',2)
plot(eta_c,tN_c,'or','linewidth',2)
l = title('Dominant mode number'); set(l,'interpreter','latex');
set(gca,'yscale','log');
l = xlabel('$\eta$'); set(l,'interpreter','latex');
l = ylabel('$\bar{t}_N$'); set(l,'interpreter','latex');
setLatexLabels2('y')
colorbar
%%
figure(23); clf
pp = 1/2/pi;

for i=1:numel(tN_c)
    tN = tN_c(i);
    eta = eta_c(i);

    Ka = 1;
    Zc  =       1;
    Ga   =      s;
    Z_Lbar = (1+eta)/(1-eta);
    delay = (1-exp(-tN*s));
    ga_i = (1+Z_Lbar*tanh(Ga))./(Z_Lbar+1*tanh(Ga));
    Ga_i = ga_i .*Ka./s.*delay;
    %%

    subplot(211);
    plot(w*pp,db([Ga_i]));
    hold on;


    subplot(212);
    plot(w*pp,radtodeg(phase(Ga_i)) + ... 
        (w>1/pp/tN)*360*0 + (w>2/pp/tN)*360*1 + (w>3/pp/tN)*360*1  );
    hold on

end

subplot(211);
plot(w*pp,w*db(1),':k');
xlim([0 max(w*pp)]);
ylabel('Magnitude (dB)');
ylim([-60 20])

subplot(212);
plot(w*pp,-180+0*w,':k'); 
ylabel('Phase (deg)');
xlabel('Frequency (Hz)');
% setLatexLabels;
xlim([0 max(w*pp)])
% l = title('$t_N=0.15$'); set(l,'interpreter','latex');
% l = legend('$\eta=0.6$','$\eta=0.2$'); set(l,'interpreter','latex');
