function be = eval_elem(a,elem,node,V1,V2,varargin)
% function be = eval_elem(a,elem,node,V1,V2,varargin)

xnode = node(elem);
gauss = calc_gauss(a,elem,3);

xgauss = gauss.coord;
N = calc_N(elem,xnode,xgauss);
[DN,detJ] = calc_DN(elem,xnode,xgauss);

c = mylocnode(a.c,elem,node,N);
tau = mylocnode(a.tau,elem,node,N);
noempty = ~isempty(V1) && ~isempty(V2);

if isempty(V1)
    V1 = N;
else
    V1 = N*localize(elem,V1);
end

if isempty(V2)
    V2 = c'*DN;
else
    V2 = c'*DN*localize(elem,V2);
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
    if size(c,1)==1
        c = c';
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

