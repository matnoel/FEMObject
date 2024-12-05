function be = eval_elem(a,elem,node,V1,V2,varargin)
% function be = eval_elem(a,elem,node,V1,V2,varargin)

xnode = node(elem);
gauss = calc_gauss(a,elem,3);

xgauss = gauss.coord;
N  = calc_N(elem,xnode,xgauss);
[DN,detJ] = calc_DN(elem,xnode,xgauss);

c = mylocnode(a.c,elem,node,N);
if a.type==1
    k = mylocnode(a.k,elem,node,N);
end

noempty = ~isempty(V1) && ~isempty(V2);

if a.type >=1
    ctemp = c;
else
    ctemp = 1;
end

if isempty(V1)
    V1 = ctemp'*DN;
else
    V1 = ctemp'*DN*localize(elem,V1);
end
if isempty(V2)
    V2 = ctemp'*DN;
else
    V2 = ctemp'*DN*localize(elem,V2);
end

switch a.type
    case 0
        cxi = norm(c);
        he = calc_h(elem,node);
        tau = he*cxi/2;
    case 1
        cxi = norm(c);
        he = 2/sum(abs(c'*DN/cxi));
        Pe = (cxi*he/2./k);
        tau = he/2./cxi.*(1./tanh(Pe)-1./Pe);
    case 2
        cxi = norm(c);
        he = 2/sum(abs(c'*DN/cxi));
        tau = he/2/cxi;
    otherwise
        error('pas prevu')
end


be = tau*V1'*V2;

be = sum(gauss.w*abs(detJ)*be,4);

if noempty
    be = double(sum(be,3));
end

