function be = vector(a,elem,node,varargin)
% function be = vector(a,elem,node,varargin)


xnode = node(elem);
gauss = calc_gauss(a,elem);
xgauss = gauss.coord;

if any([a.q,a.pk]==0)
N  = calc_N(elem,xnode,xgauss);
end
if any([a.q,a.pk]==1)
[DN,detJ] = calc_DN(elem,xnode,xgauss);
else
detJ = calc_detJ(elem,xnode,xgauss);   
end

if isempty(a.pk)
        k = a.k;
elseif a.pk==0
        k = N*localize(elem,a.k);
elseif a.pk==1
k = DN*localize(elem,a.k);
end

switch a.q
    case 0
    be = N'*k;
    case 1
    be = DN'*k;
    otherwise
        error('pas defini')
end

be = sum(gauss.w*abs(detJ)*be,4);

