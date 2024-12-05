function z=end(u,k,n)

switch n
    case 1
z=prod(size(u.X));
    case 2
z=size(u.X,k);        
end