function gauss=calc_gausspoints(h,n)
% function gauss=calc_gausspoints(h,n)
% calcul des points de gauss par morceau pour une base FE sur [0,1]
%

param = get(h,'param');
ns = 2^(param.n+1);
% p = param.p;
gausssub = calc_gausspoints(POLYLEGENDRE(),n);

gauss.coord = zeros(n,ns);
gauss.w = zeros(n,ns);

for i=1:ns
    x1 = 2^-(param.n+1)*(i-1);
    x2 = 2^-(param.n+1)*(i);
    dx = 2^-(param.n+1);
    gauss.coord(:,i) = (x1+x2)/2+gausssub.coord*dx/2;
    gauss.w(:,i) = dx*gausssub.w;
end

gauss.coord = gauss.coord(:)';
gauss.w = gauss.w(:)';
gauss.nbgauss = length(gauss.w);
