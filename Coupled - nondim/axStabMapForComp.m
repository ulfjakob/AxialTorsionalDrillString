function axStabMapForComp(p)

% nw = 1e3;
nw = 3e3;
w = linspace(3e-3,10*2*pi,nw)';
% Angular frequency the Laplace variable is evaluated at
% w = linspace(3e-3,15*2*pi,nw)';
s = 1j*w;
%%

V0Vec   = logspace(-1,1.5,60);
OmVec   = logspace(-1.31,1.71,60);

%%
Z = nan(numel(V0Vec),numel(OmVec));
Z_wu = nan(numel(V0Vec),numel(OmVec));
for i = 1:numel(V0Vec)
    i
    V0 = V0Vec(i);
    for j = 1:numel(OmVec)
        Om = OmVec(j);
        tN = 1/Om;
        Z_Lbar = (1+p.eta_a)/(1-p.eta_a);
        
        Zc = 1;
        Ga = s/p.cBar;
        delay = (1-exp(-tN*s));
        ga_i = (1+Z_Lbar*tanh(Ga))./(Z_Lbar+1*tanh(Ga));
        Ga_i = ga_i .*1./s.* p.K_a .* delay;
        
        Zc = 1;
        Ga = s;
        gt_i = (1+Z_Lbar*tanh(Ga))./(Z_Lbar+1*tanh(Ga));
        Gt_i = -gt_i .*1./s.* V0/Om .* delay;
        
        G = Ga_i + Gt_i;

        [unstab,dG,wu] = nyqStab(G,w);
        
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
figure(1); clf;
% contourf(OmVec,V0Vec,abs(Z),[0 1 2 4 8 16],'ShowText','on')
contourf(OmVec,V0Vec,db(Z+1e-2),[-40 -20:5:60],'ShowText','on')

tString = sprintf('$K_a=%d$',p.K_a);
l = title(tString); set(l,'interpreter','latex');

set(l,'interpreter','latex');
set(gca,'xscale','log')
set(gca,'yscale','log')
l = ylabel('$V_0$'); set(l,'interpreter','latex');
l = xlabel('$\Omega_0=1/t_N$'); set(l,'interpreter','latex');
xlim([0.05 40.45])
% setLatexLabels2()
hold on