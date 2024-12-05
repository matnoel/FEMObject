function b = calc_vector(a,S,varargin)
% function b = calc_matrix(a,S,varargin)
% a : LINFORM

be = cell(1,getnbgroupelem(S));
if ~isempty(a.pk)
    a.k = unfreevector(S,a.k);
end

for p=1:getnbgroupelem(S)
    be{p} = vector(a,getgroupelem(S,p),getnode(S),varargin{:});
end
b = assemble_vectorelem(S,be);

b = a.fact*b;

b = freevector(S,b);
