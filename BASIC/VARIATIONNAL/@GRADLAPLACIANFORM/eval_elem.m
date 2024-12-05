function be = eval_elem(a,elem,node,V1,V2,varargin)
% function be = eval_elem(a,elem,node,V1,V2,varargin)

xnode = node(elem);
gauss = calc_gauss(a,elem);
xgauss = gauss.coord;
N = calc_N(elem,xnode,xgauss);
[DN,detJ] = calc_DN(elem,xnode,xgauss);
DeltaN = calc_DeltaN(elem,xnode,xgauss);

c1 = mylocnode(a.c1,elem,node,N);
c2 = mylocnode(a.c2,elem,node,N);

tau = mylocelem(a.tau,elem);
noempty = ~isempty(V1) && ~isempty(V2);

if isempty(V1)
    V1 = c1'*DN;
else
    V1 = c1'*DN*localize(elem,V1);
end
if isempty(V2)
    V2 = c2'*DeltaN;
else
    V2 = c2'*DeltaN*localize(elem,V2);
end

be = tau*V1'*V2;

be = sum(gauss.w*abs(detJ)*be,4);

if noempty
    be = double(sum(be,3));
end

function c = mylocnode(c,elem,node,N)
% function c = mylocnode(c,elem,node,N)

c = full(c);
if isempty(c)
    c = 1;
elseif numel(c)~=getdim(elem)
    if size(c,1)==getnbnode(node)
        c = N*localize(elem,c);
    else
        c = mylocelem(c,elem);
        if size(c,1)==getdim(elem)^2
            c = reshape(c,[getdim(elem),getdim(elem),sizeND(c)]);
        end
    end
end
return


function tau = mylocelem(tau,elem)
% function tau = mylocelem(tau,elem)

if numel(tau)==1
    tau = MYDOUBLEND(tau);
else
    tau = tau(getnumber(elem),:);
    tau = reshape(full(tau)',size(tau,2),1,getnbelem(elem));
    tau = MYDOUBLEND(tau);
end

return

