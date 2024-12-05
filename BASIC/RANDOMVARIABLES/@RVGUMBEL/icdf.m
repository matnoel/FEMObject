function P = icdf(u,x)

param=struct2cell(getparam(u));
P = gumbelinv(x,param{:});
