function an = mynorm(a,varargin)

an = sqrt(abs(prodscal(a,a,varargin{:})));
