function P = pdf(u,x)

param=struct2cell(getparam(u));
P = gumbelpdf(x,param{:});
