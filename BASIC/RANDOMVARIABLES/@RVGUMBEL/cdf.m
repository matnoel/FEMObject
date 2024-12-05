function P = cdf(u,x)

param=struct2cell(getparam(u));
P = gumbelcdf(x,param{:});
