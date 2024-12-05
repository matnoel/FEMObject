function N = getN(elem,xi)
% function N = getN(elem,xi)

if nargin==2
    xi = double(xi);
    N = zeros([size(xi,1),3,sizeND(xi)]);
    N(:,1,:) = 1/2*(xi(:,1,:)-1).*(xi(:,1,:));
    N(:,2,:) = 1/2*(1+xi(:,1,:)).*(xi(:,1,:));
    N(:,3,:) = (1-xi(:,1,:)).*(1+xi(:,1,:));
else
    % N = inline('1/2*[(xi(1)-1).*(xi(1)),(1+xi(1)).*(xi(1)),2*(1+xi(1)).*(1-xi(1))]','xi');
    N = @(xi) 1/2*[(xi(1)-1).*(xi(1)),(1+xi(1)).*(xi(1)),2*(1+xi(1)).*(1-xi(1))];
end
