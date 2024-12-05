function apc = project(a,pc,varargin)
% function apc = project(a,pc,varargin)

pc = getPC(pc);
rv = RANDVARS(a);
[ok,rep]=ismember(rv,RANDVARS(pc));
if ~all(ok)
    error('le chaos ne correspond pas avec les variables aleatoires de la RANDVARFUNTION')
end

pca = restrictdim(pc,rep);

apc = decompfun(pca,[],[],@(x) randomeval(a,x,RANDVARS(pca)));

if getM(rv)<getM(pc)
    apc = project(apc,pc);
end
