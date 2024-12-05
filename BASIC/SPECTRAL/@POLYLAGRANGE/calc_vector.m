function b = calc_vector(l,a,f,m)
% function b = calc_vector(l,a,f)
% calcul du vecteur avec
% bi = int(li(x)^a*f(x))
% function b = calc_matrix(l,a,f,m)
% m : nombre de points pour la quadrature de
% Gauss lobatto (par defaut celle associ�e a l)

n = getnbpoints(l);
if nargin<4
    m = n+1;
end
h = getparam(l,'orthopoly');
domain = getdomain(l);
g = calc_gausslobattopoints(h,m,domain);

if a==1
    funi = @dpolyval;
elseif a==0
    funi = @polyval;
end

b = zeros(n,1);
for i=1:n
    b(i) = g.w*(funi(l,i-1,g.coord).*f(g.coord));
end

