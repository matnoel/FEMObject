function gauss=calc_gausslobattopoints(h,n,domain)

% function gauss=calc_gausslobattopoints(h,n)
%
% calcul des points de gauss associes a la mesure uniforme 1/2 sur [-1,1]
% n : nombre de points (incluant les extrémités)
% 
% function gauss=calc_gausslobattopoints(h,n,domain)
% calcul des points de gauss associes a la mesure uniforme 1 sur [domain(1),domain(2)]

if n<2
    error('2 points de Gauss')
end

m = n-1;

gauss.w=zeros(1,n);
gauss.coord=zeros(n,1);

hprimem = -[0,m*polycoeff(h,m)/sqrt(2*m+1)]+[m*polycoeff(h,m-1)/sqrt(2*m-1),0,0]; 

x=roots(fliplr(hprimem));
x=sort(x);
gauss.coord=x(:) ;

l = POLYLAGRANGE(x);
g = calc_gausspoints(h,n);
g.w = g.w;

for i=1:n
gauss.w(i) = g.w(:)'*polyval(l,i-1,g.coord(:));
end

gauss.nbgauss=n;

if nargin==3    
   dx = domain(2)-domain(1);
   gauss.coord = transfer(RVUNIFORM(-1,1),RVUNIFORM(domain(1),domain(2)),gauss.coord);
   gauss.w = gauss.w*dx;
end
