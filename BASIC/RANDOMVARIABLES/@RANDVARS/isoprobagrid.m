function A = isoprobagrid(u,P,varargin)

r = RANDVARS(RVUNIFORM(0,1),u.M);

if P>1
n = P-1;
else    
p = P^(1/u.M);
n = ceil(1/p);
end

for k =1 : u.M
    b{k} = linspace(1/n,1-1/n,n);
end

A = cell(1,u.M);
[A{:}] = ndgrid(b{:});

for k =1 : u.M
   A{k}=A{k}(:); 
end

A = transfer(r,u,A);
