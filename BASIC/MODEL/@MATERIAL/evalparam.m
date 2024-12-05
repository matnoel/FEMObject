function varargout = evalparam(mat,name,varargin)
% function [val,x] = evalparam(mat,name,elem,xnode,xgauss)

param = getparam(mat,name);

if nargin==3
    x = varargin{1};
else
    elem = varargin{1};
    xnode = varargin{2};
    xgauss = varargin{3};
    if isa(param,'inline') || isa(param,'function_handle') || isa(param,'FENODEFIELD')
        [x,Nuni] = calc_x(elem,xnode,xgauss);
    end
end

if isa(param,'FENODEFIELD')
    % param = full(double(param));
    % ss = size(getconnec(elem)');
    % param = param(getconnec(elem)',:);
    % 
    % param = reshape(param,[ss(1) ss(2) size(param,2)]);
    % keyboard
    % param = permute(param,[3,1,2]);
    % val = Nuni*MYDOUBLEND(param);
    % 
    param = full(double(param));
    param = param(getconnec(elem)');
    param = reshape(param,[size(param,1) 1 size(param,2)]);
    
    val = Nuni*MYDOUBLEND(param);

elseif isa(param,'double') || isa(param,'MYDOUBLE')
    val = double(param);
    % sx = size3D(x);
    % val = repmat(param,[sx(1) 1  sx(3)]);
    
elseif isa(param,'MYDOUBLEND')
    val = param;
    
elseif isa(param,'FEELEMFIELD')
    val = MYDOUBLEND(param);
    
elseif isa(param,'inline') || isa(param,'function_handle')
    sx = size(x);
    ndim = [1,3:length(sx),2];
    sxbis = sx(ndim);
    param = fcnchk(param); % transformation de la fonction en inline (fonction matlab)
    val = feval(param,reshape(permute(x,ndim),prod(sxbis(1:end-1)),sxbis(end)));
    val = MYDOUBLEND(reshape(val,[sxbis(1) 1 sxbis(2:end-1)]));
else
    error('param non defini')
end

varargout{1} = val;
if nargout>1
    varargout{2} = x;
end
