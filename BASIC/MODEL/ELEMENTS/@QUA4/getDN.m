function DN = getDN(elem,xi)
% function DN = getDN(elem,xi)

if nargin==2
    xi = double(xi);
    DN=zeros([2,4,sizeND(xi)]);
    eta = xi(:,2,:);
    xi = xi(:,1,:);
    
    DN(1,1,:) = -1/4*(1-eta);
    DN(1,2,:) = 1/4*(1-eta);
    DN(1,3,:) = 1/4*(1+eta);
    DN(1,4,:) = -1/4*(1+eta);
    DN(2,1,:) = -1/4*(1-xi);
    DN(2,2,:) = -1/4*(1+xi);
    DN(2,3,:) = 1/4*(1+xi);
    DN(2,4,:) = 1/4*(1-xi);
else
    % DN = inline('1/4*[-(1-xi(2)), (1-xi(2)), (1+xi(2)),-(1+xi(2)); -(1-xi(1)), -(1+xi(1)), (1+xi(1)),(1-xi(1))]','xi');
    DN = @(xi) 1/4*[-(1-xi(2)), (1-xi(2)), (1+xi(2)),-(1+xi(2)); -(1-xi(1)), -(1+xi(1)), (1+xi(1)),(1-xi(1))];
end
