%%
clc
clear;
% close all;

% nw = 1e3;
nw = 3e3;
% nw = 1e4;
w = linspace(3e-3,10*2*pi,nw)';
% Angular frequency the Laplace variable is evaluated at
% w = linspace(3e-3,15*2*pi,nw)';
s = 1j*w;

%% Physical Params
eta     = .7;
V0Vect  = logspace(-1,1.5,120);
OmVec   = logspace(-1.31,1.61,120);
% V0Vec   = logspace(-1,1.5,60);
% OmVec   = logspace(-1.31,1.61,60);

%%
Z = nan(numel(OmVec),numel(V0Vect));
Z_wu = nan(numel(OmVec),numel(V0Vect));
for i = 1:numel(V0Vect)
    i
    V0 = V0Vect(i);
    for j = 1:numel(OmVec)
        Om = OmVec(j);
        tN = 1/Om;
        Z_Lbar = (1+eta)/(1-eta);
        
        Zc = 1;
        Ga = s;
        delay = (1-exp(-tN*s));
        gt_i = (1+Z_Lbar*tanh(Ga))./(Z_Lbar+1*tanh(Ga));
        Gt_i = -gt_i .*1./s.* V0/Om .* delay;

        [unstab,dG,wu] = nyqStab(Gt_i,w);
        
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
figure(2); clf;
% contourf(OmVec,V0Vec,db(Z),'ShowText','on')
% contourf(OmVec,V0Vec,db(Z),[-70 -40 -20 0 20 40 70],'ShowText','on')
Zss = Z>1;
contourf(OmVec,V0Vect,Zss,[0 .1:.2:.9 .99],'ShowText','off')
colorbar
% tn = sprintf('Inverse gain margin (dB) $M_G = \\max_{\\varpi\\in \\varpi_{180}} \\Big|G_t(j\\varpi) \\Big|, \\quad \\eta_t= %.2f$',eta);
tn = sprintf('Torsional stability, $M_T$, $\\quad \\eta_t= %.2f$',eta);
% l = title('Inverse gain margin (dB) $M_G = \max_{\varpi\in \varpi_{180}} \Big|G_t(j\varpi) \Big|, \quad \eta$');
% l = title(tn);
% set(l,'interpreter','latex');
set(gca,'xscale','log')
set(gca,'yscale','log')
l = ylabel('$V_0$'); set(l,'interpreter','latex');
l = xlabel('$\Omega_0=1/t_N$'); set(l,'interpreter','latex');
xlim([0.05 40.45])
setLatexLabels2()

%% Draw line
OmLine = nan(numel(V0Vect),1);
for i = 1:numel(V0Vect)
    OmInd = min(max(find(Zss(i,:)<0.5,1)-1,0),numel(V0Vect)-1);
%     OmLine(i) = mean(OmVec( OmInd:OmInd+1 ) );
    OmLine(i) = OmVec( OmInd );

end
hold on;

%
plot(smooth(OmLine,1),V0Vect,'-r','linewidth',3)

%%
% h = figure(2);
% if eta < 0
% %     tS = sprintf('b)');
% %     l = title(tS); set(l,'interpreter','latex','fontsize',20,...
% %     'Units', 'normalized','Position', [0.01 1], 'HorizontalAlignment', 'left');
%     pdfmatlabfrag2(h,'plotSaves/stabMap_NEG');
% else
% %     tS = sprintf('a)');
% %     l = title(tS); set(l,'interpreter','latex','fontsize',20,...
% %     'Units', 'normalized','Position', [0.01 1], 'HorizontalAlignment', 'left');
%     pdfmatlabfrag2(h,'plotSaves/stabMap');
% end















