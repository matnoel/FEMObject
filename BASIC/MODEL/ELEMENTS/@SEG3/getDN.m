function DN = getDN(elem,xi)
% function DN = getDN(elem,xi)

if nargin==2
    xi = double(xi);
    DN = zeros([1,3,sizeND(xi)]);
    xi = xi(:,1,:);
    DN(1,1,:) = xi-1/2;
    DN(1,2,:) = 1/2+xi;
    DN(1,3,:) = -2*xi;
end
