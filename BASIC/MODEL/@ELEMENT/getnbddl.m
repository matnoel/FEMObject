function n = getnbddl(elem,varargin)
% function n = getnbddl(elem,varargin)

if nargin==1
    n = elem.nbddl;
else
    np = getnbddlpernode(elem,varargin{:});
    n  = getnbnode(elem)*np;
end


