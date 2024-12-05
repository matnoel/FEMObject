function Ae = matrix(a,elem,node,varargin)
% function Ae = matrix(a,elem,node,varargin)

xnode = node(elem);
gauss = calc_gauss(a,elem);
xgauss = gauss.coord;

if any([a.p,a.q,a.pk]==0)
N  = calc_N(elem,xnode,xgauss);
end
if any([a.p,a.q,a.pk]==1)
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


if a.p==1 && a.q==1
    Ae = k.*(DN'*DN);
elseif a.p==0 && a.q==0
    Ae = k.*N'*N;        
elseif a.p==1 && a.q==0
    Ae = N'*(k'*DN);             
elseif a.p==0 && a.q==1
    Ae = (k'*DN)'*N;
end

Ae = sum(gauss.w*abs(detJ)*Ae,4);

