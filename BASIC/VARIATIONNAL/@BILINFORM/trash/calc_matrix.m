function A = calc_matrix(a,S,varargin)
% function A = calc_matrix(a,varargin)
% a : BILINFORM

Ae = cell(1,getnbgroupelem(S));
if ~isempty(a.pk)
    a.k = unfreevector(S,a.k);
end

for p=1:getnbgroupelem(S)
    Ae{p} = matrix(a,getgroupelem(S,p),getnode(S),varargin{:});
end
A = assemble_matrixelem(S,Ae);

A = freematrix(S,A);

A = a.fact*A;
