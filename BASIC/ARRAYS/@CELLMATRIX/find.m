function [i,j,k,s] = find(u)
% function [i,j,k,s] = find(u)

switch nargout
    case {0,1}
        i = find(u.value) ;
    case 2
        [i,j] = find(u.value) ;
    case {3,4}
        [i,k,s] = find(u.value) ;
        [i,j] = ind2sub(size(u),i) ;
end