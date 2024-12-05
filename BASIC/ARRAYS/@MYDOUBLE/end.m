function z = end(u,k,n)
% function z = end(u,k,n)

switch n
    case 1
        z = prod(size(u));
    otherwise
        z = size(u,k);
end
