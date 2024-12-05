function b = calc_vector2D(S,a,f,mx,my)
% function b = calc_vector2D(S,a,f)
% calcul du vecteur avec
% bi = int(li(x)^a*f(x))
% function b = calc_matrix(l,a,f,m)
% m : nombre de points pour la quadrature de 
% Gauss lobatto (par defaut celle associée a l)

lx = S.L;
ly = S.L;
nx = getnbpoints(lx);
ny = getnbpoints(ly);
N = nx*ny;
[I,J]=meshgrid(1:nx,1:ny);

if nargin<6
mx = getnbpointsperelement(lx)+1;
my = getnbpointsperelement(ly)+1;
end
gx = calc_gausslobattopoints(lx,mx);
gy = calc_gausslobattopoints(ly,my);
g = tensorize_quadrature_rule({gx,gy});


if a==0
b = zeros(N,1);
for i=1:N
    lval = polyval(lx,I(i)-1,g.coord(:,1)).*polyval(ly,J(i)-1,g.coord(:,2));
    b(i)=g.w*(lval.*f(g.coord));
end
else
    error('pas programme')
end

