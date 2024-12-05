function u = orth(u,dim)
% function u = orth(u,dim)
% apply the function orth to all cells u.F{i,j} with j in dim
% if nargin=1, all dimensions j are considered

if nargin==1
    dim = 1:u.dim;
end

u.F(:,dim) = cellfun(@(C) orth(full(C)),u.F(:,dim),'UniformOutput',false);
% for i=1:u.m
%  for k=dim
%    u.F{i,k}=orth(full(u.F{i,k})); 
%  end
% end

