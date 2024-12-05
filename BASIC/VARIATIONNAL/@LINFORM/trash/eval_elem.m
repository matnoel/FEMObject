function be = eval_elem(a,V,elem,node,varargin)
% function be = eval_elem(a,V,elem,node,varargin)

xnode = node(elem);
gauss = calc_gauss(a,elem);
xgauss = gauss.coord;

if any([a.q,a.pk]==0)
    N = calc_N(elem,xnode,xgauss);
end
if any([a.q,a.pk]==1)
    [DN,detJ] = calc_DN(elem,xnode,xgauss);
else
    detJ = calc_detJ(elem,xnode,xgauss);
end

if isempty(a.pk)
    k = a.k;
else
    switch a.pk
        case 0
            k = N*localize(elem,a.k);
        case 1
            k = DN*localize(elem,a.k);
    end
end

Ve = localize(V,elem);

if a.q==1
    DV = DN*Ve;
    be = DV'*k;
elseif a.q==0
    V = N*Ve;
    be = V.*k;
end

be = sum(gauss.w*abs(detJ)*be,4);
be = sum(be,3);



