%% MGVec = gainMargPipe(tN_vec,eta_vec,Ka_vec,etaType)
%
% Returns a vector of gain margins for the parameter vectors parsed as
% arguments.
%
function MGVec = gainMargPipe(tN_vec,eta_vec,Ka_vec,etaType)

if nargin == 3
    etaType = 'k';
end

% Angular frequency the Laplace variable is evaluated at
nw = 1e4;
w = linspace(3e-3,10*4*pi,nw)';
s = 1j*w;
cBar = 1.6;
%%
MGVec = nan(numel(tN_vec),1);
for i=1:numel(tN_vec)
    tN = tN_vec(i);
    eta = eta_vec(i);
    Ka = Ka_vec(i);
    
    if strcmp(etaType,'k')
        k = -log(abs(eta));
        Zc  =      sqrt(1+k./s);
        Ga  =      s/cBar.*sqrt(1+k./s);
        if eta > 0
            ga_i = 1./Zc.*tanh(Ga);
        else
            ga_i = 1./Zc./tanh(Ga);
        end
    elseif strcmp(etaType,'Z')
        Zc  =      1;
        Ga  =      s/cBar;
        Z_Lbar = (1+eta)/(1-eta);
        ga_i = (1+Z_Lbar*tanh(Ga))./(Z_Lbar+1*tanh(Ga));
    end
    
    delay = (1-exp(-tN*s));
    Ga_i = ga_i .*Ka./s.*delay;

    [~,MGVec(i),~] = nyqStab(Ga_i,w);

end