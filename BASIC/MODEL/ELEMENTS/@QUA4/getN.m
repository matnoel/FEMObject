function N = getN(elem,xi)
% function N = getN(elem,xi)

if nargin==2
    xi = double(xi);
    N = zeros([size(xi,1),4,sizeND(xi)]);
    N(:,1,:) = 1/4*(1-xi(:,1,:)).*(1-xi(:,2,:));
    N(:,2,:) = 1/4*(1+xi(:,1,:)).*(1-xi(:,2,:));
    N(:,3,:) = 1/4*(1+xi(:,1,:)).*(1+xi(:,2,:));
    N(:,4,:) = 1/4*(1-xi(:,1,:)).*(1+xi(:,2,:));
else
    % N = inline('1/4*[(1-xi(1)).*(1-xi(2)),(1+xi(1)).*(1-xi(2)),(1+xi(1)).*(1+xi(2)),(1-xi(1)).*(1+xi(2))]','xi');
    N = @(xi) 1/4*[(1-xi(1)).*(1-xi(2)),(1+xi(1)).*(1-xi(2)),(1+xi(1)).*(1+xi(2)),(1-xi(1)).*(1+xi(2))];
end
