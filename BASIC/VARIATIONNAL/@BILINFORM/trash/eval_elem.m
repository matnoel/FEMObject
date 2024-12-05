function Ae = eval_elem(a,U,V,elem,node,varargin)
% function Ae = eval_elem(a,U,V,elem,node,varargin)

xnode = node(elem);
gauss = calc_gauss(a,elem);
xgauss = gauss.coord;

if any([a.p,a.q,a.pk]==0)
    N = calc_N(elem,xnode,xgauss);
end
if any([a.p,a.q,a.pk]==1)
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

Ue = localize(U,elem);
Ve = localize(V,elem);

if a.p==1 && a.q==1
    DU = DN*Ue;
    DV = DN*Ve;
    Ae = k.*(DV'*DU);
elseif a.p==0 && a.q==0
    U = N*Ue;
    V = N*Ve;
    Ae = k.*V'*U;
elseif a.p==1 && a.q==0
    DU = DN*Ue;
    V = N*Ve;
    Ae = V'*(k'*DU);
elseif a.p==0 && a.q==1
    U = N*Ue;
    DV = DN*Ve;
    Ae = (k'*DV)'*U;
end

Ae = sum(gauss.w*abs(detJ)*Ae,4);
Ae = sum(Ae,3);



