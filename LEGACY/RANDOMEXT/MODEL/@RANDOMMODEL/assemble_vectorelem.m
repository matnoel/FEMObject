function f = assemble_vectorelem(S,fe,varargin)
% function f = assemble_vectorelem(S,fe,varargin)

if S.nbgroupelem==0
    f = sparse(S.nbddl,1);
else
    f = sparse(S.nbddl,size(fe{1},2));
end

for p=1:S.nbgroupelem
    [ip,valp] = vectorelem(S.groupelem{p},double(fe{p}));
    f(ip,:) = f(ip,:)+valp;
end



