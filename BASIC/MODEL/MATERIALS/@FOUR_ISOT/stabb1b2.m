function ke = stabb1b2(mat,elem,xnode,xgauss,varargin)
% function ke = stabb1b2(mat,elem,xnode,xgauss,varargin)

k = evalparam(mat,'k',elem,xnode,xgauss);
B = calc_B(elem,xnode,xgauss);
b = evalparam(mat,'b',elem,xnode,xgauss);
bstab = evalparam(mat,'bstab',elem,xnode,xgauss);

bxi = norm(bstab);
he = 2/sum(abs(bstab'*B/bxi));
if getparam(mat,'stabilize')==1
    Pe = (bxi*he/2/k);
    tau = he/2./bxi.*(1./tanh(Pe)-1./Pe);
elseif getparam(mat,'stabilize')==2
    tau = he/2/bxi;
end
ke = tau*(B'*bbstab)*(b'*B);
