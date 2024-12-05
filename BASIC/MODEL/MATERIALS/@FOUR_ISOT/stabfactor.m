function tau = stabfactor(mat,elem,xnode,xgauss,varargin)
% function tau = stabfactor(mat,elem,xnode,xgauss,varargin)

if getparam(mat,'stabilize')==0
    mat = setparam(mat,'stabilize',1);
end
k = evalparam(mat,'kstab',elem,xnode,xgauss);
B = calc_B(elem,xnode,xgauss);
bstab = evalparam(mat,'bstab',elem,xnode,xgauss);

bxi = norm(bstab);
he = 2/sum(abs(bstab'*B/bxi));
if getparam(mat,'stabilize')==1
    Pe = (bxi*he/2/k);
    tau = he/2./bxi.*(1./tanh(Pe)-1./Pe);
elseif getparam(mat,'stabilize')==2
    tau = he/2/bxi;
end
