function DDN = getDDN(elem,xi)
% function DDN = getDDN(elem,xi)

if nargin==2
    DDN = repmat(zeros(3,3),[1,1,sizeND(xi)]);    
else
    % DDN = inline('repmat(zeros(3,3),[1,1,sizeND(xi)])','xi');
    DDN = @(xi) repmat(zeros(3,3),[1,1,sizeND(xi)]);
end
