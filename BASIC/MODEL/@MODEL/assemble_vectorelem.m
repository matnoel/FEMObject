function f = assemble_vectorelem(S,fe,varargin)
% function f = assemble_vectorelem(S,fe,varargin)

liste = getcharin('selgroup',varargin,1:S.nbgroupelem);
if S.nbgroupelem==0
    f = sparse(S.nbddl,1);
else
    f = sparse(S.nbddl,size(fe{liste(1)},2));
end

for p=liste
    [ip,valp] = vectorelem(S.groupelem{p},double(fe{p}));
    f(ip,:) = f(ip,:)+valp;
end



