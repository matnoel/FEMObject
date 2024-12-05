function apc = project(a,pc,varargin)
% function apc = project(a,pc,varargin)

stodim = getstodim(a);
[ok,rep]=ismember([stodim{:}],pc);
if ~all(ok)
    error('la variable conditionnelle n''est pas incluse dans le chaos')
elseif length(stodim)~=length([stodim{:}])
    error('toutes les variables doivent etre numerotees')
end

PC = restrictdim(pc,[stodim{:}]);

ng = getcharin('nbgauss',varargin);
apc = decompfun(PC,ng,[],@(x,RV,a) condtransfer(RV,a,x),RANDVARS(PC),a);

apc = project(apc,pc);
