function ke = stab(mat,elem,xnode,xgauss,varargin)
% function ke = stab(mat,elem,xnode,xgauss,varargin)

k = evalparam(mat,'D2',elem,xnode,xgauss);
B = calc_B(elem,xnode,xgauss);
b = evalparam(mat,'D1',elem,xnode,xgauss);

if getparam(mat,'stabilize')>0
    bxi = norm(b);
    he = 2/sum(abs(b'*B/bxi));
    if getparam(mat,'stabilize')==1
        Pe = (bxi*he/2/k);
        tau = he/2./bxi.*(1./tanh(Pe)-1./Pe);
    elseif getparam(mat,'stabilize')==2
        tau = he/2/bxi;
    end
    ke = tau*(B'*b)*(b'*B);
    
end
