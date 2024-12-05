function ke = reac(mat,elem,xnode,xgauss,varargin)
% function ke = reac(mat,elem,xnode,xgauss,varargin)

if nargin>4 && isa(varargin{1},'POLYCHAOS')
    r = evalparampc(mat,'r',varargin{1},elem,xnode,xgauss);
else
    r = evalparam(mat,'r',elem,xnode,xgauss);
end
N = calc_N(elem,xnode,xgauss);
ke = N'*full(r)*N;
