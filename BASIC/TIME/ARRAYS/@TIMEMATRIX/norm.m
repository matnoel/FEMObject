function an = norm(a,varargin)
% function an = norm(a,varargin)

an = sqrt(full(abs(prodscal(a,a,varargin{:}))));
