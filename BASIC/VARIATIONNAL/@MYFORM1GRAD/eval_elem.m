function be = eval_elem(a,elem,node,V,W,U1,U2,varargin)
% function be = eval_elem(a,elem,node,V,W,U1,U2,varargin)

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



U1loc = localize(elem,U1);
U2loc = localize(elem,U2);
DU1 = DN*U1loc;
DU2 = DN*U2loc;
if ~isempty(W)
    Wloc = localize(elem,W);
    DW = DN*Wloc;
end

if approx
    U1U2 = N*localize(elem,U1.*U2);
    U1N = U1loc'.*N;
    U2N = U2loc'.*N;
    if ~isempty(W)
        U1W = N*localize(elem,U1.*W);
        U2W = N*localize(elem,U2.*W);
    end
else
    U1 = N*U1loc;
    U2 = N*U2loc;
    U1U2 = U1.*U2;
    if ~isempty(W)
        W = N*Wloc;
        U1W = U1.*W;
        U2W = U2.*W;
    end
end

if isempty(V) && isempty(W)
    if approx==0
        be = DN'*(DN*U1U2) + DN'*(DU1*U2 + DU2*U1)*N ;
    else
        be = DN'*(DN*U1U2) + DN'*(DU1*U2N + DU2*U1N) ;
    end
    be = sum(gauss.w*abs(detJ)*be,4);
elseif isempty(V)
    be = DN'*(DW*U1U2) + DN'*(DU1*U2W + DU2*U1W) ;
    be = sum(gauss.w*abs(detJ)*be,4);
else
    DV = DN*localize(elem,V);
    be = DV'*(DW*U1U2 + DU1*U2W + DU2*U1W) ;
    be = sum(gauss.w*abs(detJ)*be,4);
    be = double(sum(be,3));
end


