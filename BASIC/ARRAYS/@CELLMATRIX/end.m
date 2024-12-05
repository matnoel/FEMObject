function z=end(u,k,n)

switch n
    case 1
z=prod(u.s);
    case 2
z=u.s(k);        
end