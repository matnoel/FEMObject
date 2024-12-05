function DDN = getDDN(elem,xi)
% function DDN = getDDN(elem,xi)

if nargin==2
    DDN = repmat([4,4,0,-8,0,0;4,0,4,0,0,-8;4,0,0,-4,4,-4],[1,1,sizeND(xi)]);
else
    % DDN = inline('repmat([4,4,0,-8,0,0;4,0,4,0,0,-8;4,0,0,-4,4,-4],[1,1,sizeND(xi)])','xi');
    DDN = @(xi) repmat([4,4,0,-8,0,0;4,0,4,0,0,-8;4,0,0,-4,4,-4],[1,1,sizeND(xi)]);
end
