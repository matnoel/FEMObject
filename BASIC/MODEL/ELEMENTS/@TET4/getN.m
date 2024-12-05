function N = getN(elem,xi)
% function N = getN(elem,xi)

if nargin==2
    ismydoublend = isa(xi,'MYDOUBLEND');
    xi = double(xi);
    N = zeros([size(xi,1),4,sizeND(xi)]);
    theta = xi(:,3,:);
    eta = xi(:,2,:);
    xi = xi(:,1,:);
    N(:,1,:) = 1-xi-eta-theta;
    N(:,2,:) = xi;
    N(:,3,:) = eta;
    N(:,4,:) = theta;
    if ismydoublend
        N = MYDOUBLEND(N);
    end
else
    % N = inline('[1-xi(1)-xi(2)-xi(3),xi(1),xi(2),xi(3)]','xi');
    N = @(xi) [1-xi(1)-xi(2)-xi(3),xi(1),xi(2),xi(3)];
end
