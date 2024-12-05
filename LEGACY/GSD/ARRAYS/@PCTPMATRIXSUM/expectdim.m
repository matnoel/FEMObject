function x = expectdim(dim,x,varargin)
% function x = expectdim(dim,x,varargin)

nodim = setdiff(1:getnbdim(x),dim);
x = expectnodim(nodim,x,varargin{:});

