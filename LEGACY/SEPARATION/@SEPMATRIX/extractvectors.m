function V = extractvectors(u,dim)

if nargin==1
    dim = 1:u.dim;
end
V = cell(1,length(dim));
for k=1:length(dim)
    V{k} = zeros(size(u.F{1,dim(k)},1),u.m);
for i=1:u.m
    V{k}(:,i) = u.F{i,dim(k)}; 
end
end

if length(V)==1
    V=V{1};
end