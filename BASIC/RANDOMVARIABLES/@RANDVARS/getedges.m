function A = getedges(u,varargin)

A = cell(1,u.M);
b = cell(1,u.M);
for k =1 : u.M
    b{k} = getedges(u.RV{k});
end

[A{:}] = ndgrid(b{:});

for k =1 : u.M
   A{k}=A{k}(:); 
end

if ~ischarin('cell',varargin)
    A=[A{:}];
end