function be = eval_elem(a,elem,node,V1,V2,varargin)
% function be = eval_elem(a,elem,node,V1,V2,varargin)

xnode = node(elem);
gauss = calc_gauss(a,elem);
xgauss = gauss.coord;
N = calc_N(elem,xnode,xgauss);
[DN,detJ] = calc_DN(elem,xnode,xgauss);

c1 = mylocnode(a.c1,elem,node,N);
c2 = mylocnode(a.c2,elem,node,N);

if size(a.tau,1)==getnbnode(node)
    tau = mylocnode(a.tau,elem,node,N);
else
    tau = mylocelem(a.tau,elem);
end
noempty = ~isempty(V1) && ~isempty(V2);

if isempty(V1)
    V1 = c1'*DN;
else
    V1 = c1'*DN*localize(elem,V1);
end
if isempty(V2)
    V2 = c2'*DN;
else
    V2 = c2'*DN*localize(elem,V2);
end

be = tau*V1'*V2;

be = sum(gauss.w*abs(detJ)*be,4);

if noempty
    be = double(sum(be,3));
end


