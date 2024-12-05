function be = eval_elem(a,elem,node,V,U,varargin)
% function be = eval_elem(a,elem,node,V,U,varargin)

approx = getapprox(a);

xnode = node(elem);
if approx
    gauss = calc_gauss(elem,1);
else
    gauss = calc_gauss(elem,2);
end
xgauss = gauss.coord;

N = calc_N(elem,xnode,xgauss);
[DN,detJ] = calc_DN(elem,xnode,xgauss);

if approx
    U2 = N*localize(elem,U.*U);
    DU = DN*localize(elem,U);
else
    Uloc = localize(elem,U);
    DU = DN*Uloc;
    U = N*Uloc;
    U2 = U.*U;
end

if isempty(V)
    if isa(U,'PCMYDOUBLEND')
        be = DN'*([DU(1).*U2;DU(2).*U2]);
    else
        be = DN'*(DU*U2);
    end
    be = sum(gauss.w*abs(detJ)*be,4);
else
    DV = DN*localize(elem,V);
    be = (DV'*DU)*U2;
    be = sum(gauss.w*abs(detJ)*be,4);
    be = double(sum(be,3));
end


