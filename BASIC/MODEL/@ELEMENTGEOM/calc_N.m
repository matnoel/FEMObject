function [N,detJ,x] = calc_N(elem,xnode,xgauss,varargin)
% function [N,detJ,x] = calc_N(elem,xnode,xgauss,varargin)

if getparam(elem,'initializeBN')
    N = getparam(elem,'N');
    if nargin==2
        detJ = getparam(elem,'detJ');
    end
    if nargin==3
        x = getparam(elem,'x');
    end
    return
end

nbelem = getnbelem(elem);

Nuni = getN(elem,xgauss);
if size(Nuni,3)==1
    Nuni = repmat(Nuni,[1,1,getnbelem(elem)]);
end

if ischarin('type',varargin)
    type = getcharin('type',varargin,'ddlnode');
    n = getnbddlpernode(elem,type);
else
    n = getcharin('nbddlpernode',varargin,getnbddlpernode(elem,'ddlnode'));
end
nbddl = n*getnbnode(elem);

if nargout>1
    detJ = calc_detJ(elem,xnode,xgauss);
end
if nargout>2
    x = calc_x(elem,xnode,xgauss);
end

N = zerosND([n,nbddl,sizeND(Nuni)]);
for i=1:n
    N(i,i:n:end) = Nuni;
end
