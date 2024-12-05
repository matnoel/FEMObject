function b = calc_vector(l,a,f,varargin)
% function b = calc_vector(l,a,f)
% calcul du vecteur avec
% bi = int(li(x)^a*f(x))
%
% function b = calc_matrix(l,a,f,m)
% N : nombre de points par element pour la quadrature de
% Gauss lobatto (par defaut celle associ�e a l)


param = getparam(l);
n = param.n;
N = getnbpoints(l);
if nargin<4
    m = getparam(l,'m');
end

if a==1
    funi = @dpolyval;
elseif a==0
    funi = @polyval;
end

b = zeros(N,1);
for k=1:n
    rep = param.reppoints{k};
    subpoints = param.subpoints{k};
    lk = POLYLAGRANGE(subpoints);
    g = calc_gausslobattopoints(POLYLEGENDRE(),param.m+1,getdomain(lk));
    for i=1:length(rep)
        b(rep(i)) = b(rep(i))+g.w*(funi(lk,i-1,g.coord).*f(g.coord));
    end
end

