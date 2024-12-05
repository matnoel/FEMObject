function DN = getDN(elem,xi)
% function DN = getDN(elem,xi)

if nargin==2
    DN = repmat([-1,1,0;-1,0,1],[1,1,sizeND(xi)]);
else
    % DN = inline('repmat([-1,1,0;-1,0,1],[1,1,sizeND(xi)])','xi');
    DN = @(xi) repmat([-1,1,0;-1,0,1],[1,1,sizeND(xi)]);
end
