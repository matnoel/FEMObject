function N = getN(elem,xi)
% function N = getN(elem,xi)

if nargin==2
    xi = double(xi);
    N = zeros([size(xi,1),2,sizeND(xi)]);
    N(:,1,:) = (1-xi(:,1,:))/2 ;
    N(:,2,:) = (1+xi(:,1,:))/2 ;
else
    % N = inline('[1-xi,1+xi]/2','xi');
    N = @(xi) [1-xi,1+xi]/2;
end
