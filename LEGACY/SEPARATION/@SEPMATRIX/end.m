function z=end(u,k,n)

switch k
    case 1
z=size(u,1);
    case 2
z=size(u,2);
    otherwise
if n==4        
z = size(u,k-2,1);
elseif n==3
z = prod(size(u,[],1));    
end
end
