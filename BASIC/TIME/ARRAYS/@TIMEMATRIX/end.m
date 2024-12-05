function z=end(u,k,n)

switch n
    case 1
z=prod(size(u));
    case 2
z=size(u,k);        
end