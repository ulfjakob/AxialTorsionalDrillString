Tbip = zeros(size(p.t));
for i=1:numel(p.t)-1
    k = p.Nt-i;
    Tbip(k) = Tbip(k+1)*(1+tLp) + y(4,k+1)*tLp;
end