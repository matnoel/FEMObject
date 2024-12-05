function N = getN(elem,xi)
% function N = getN(elem,xi)

if nargin==2
    ismydoublend = isa(xi,'MYDOUBLEND');
    xi = double(xi);
    N = zeros([size(xi,1),3,sizeND(xi)]);
    eta = xi(:,2,:);
    xi = xi(:,1,:);
    N(:,1,:) = 1-xi-eta;
    N(:,2,:) = xi;
    N(:,3,:) = eta;
    if ismydoublend
        N = MYDOUBLEND(N);
    end
else
    % N = inline('[1-xi(1)-xi(2),xi(1),xi(2)]','xi');
    N = @(xi) [1-xi(1)-xi(2),xi(1),xi(2)];
end
