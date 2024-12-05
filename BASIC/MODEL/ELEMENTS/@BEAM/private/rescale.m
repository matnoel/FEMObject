function N = rescale(N,elem)
% function N = rescale(N,elem)

nbelem = getnbelem(elem);
param = getparam(elem);
L = param.L;
if size(N,3)==1
    N = repmat(N,[1 1 nbelem])  ;
end
N = MYDOUBLEND(N) ;
s = sizeND(N);
L = repmat(L,[1,1,1,s(2:end)]);

N(1,[2,4])=N(1,[2,4])*L;
