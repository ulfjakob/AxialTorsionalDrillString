%% [MG,omega] = ustabModes(tN,eta,Ka,etaType)
%
% Inverse gain margin and frequency of all unstable modes
%
function [MG,wu] = ustabModes(tN,eta,Ka,etaType)

if nargin == 3
    etaType = 'k';
end

% Angular frequency the Laplace variable is evaluated at
nw = 1e4;
w = linspace(3e-3,10*4*pi,nw)';
s = 1j*w;
cBar = 1.6;

%%
if strcmp(etaType,'k')
    k   =  -log(abs(eta));
    Zc  =      sqrt(1+k./s);
    Ga  =   s/cBar.*sqrt(1+k./s);
    if eta > 0
        ga_i = 1./Zc.*tanh(Ga);
    else
        ga_i = 1./Zc./tanh(Ga);
    end
elseif strcmp(etaType,'Z')
    Zc  =      1;
    Ga  =      s/cBar;
    Z_Lbar = (1+eta)/(1-eta);
    ga_i = (1+Z_Lbar*tanh(Ga))./(Z_Lbar+tanh(Ga));
end
delay = (1-exp(-tN*s));
Ga_i = ga_i .*Ka./s.*delay;

N_maxModes = 10;
[unstab,MG,wu] = nyqStab(Ga_i,w,N_maxModes);

i = 0;
while i<N_maxModes && MG(i+1)>1
    i = i+1;
end
MG = MG(1:i);
wu = wu(1:i);