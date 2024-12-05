function ke = diff(mat,elem,xnode,xgauss,varargin)
% function ke = diff(mat,elem,xnode,xgauss,varargin)

if nargin>4 && isa(varargin{1},'POLYCHAOS')
    k = evalparampc(mat,'k',varargin{1},elem,xnode,xgauss);
else
    k = evalparam(mat,'k',elem,xnode,xgauss);
end
B = calc_B(elem,xnode,xgauss);
ke = B'*full(k)*B;
