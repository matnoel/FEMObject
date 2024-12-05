function N = getN(elem,xi)
% function N = getN(elem,xi)

if nargin==2
    ismydoublend = isa(xi,'MYDOUBLEND');
    xi = double(xi);
    N = zeros([size(xi,1),10,sizeND(xi)]);
    theta = xi(:,3,:);
    eta = xi(:,2,:);
    xi = xi(:,1,:);
    N(:,1,:) = 1+xi.*(2*xi+4*eta-3)+eta.*(2*eta+4*theta-3)+theta.*(2*theta+4*xi-3);
    N(:,2,:) = xi.*(2*xi-1);
    N(:,3,:) = eta.*(2*eta-1);
    N(:,4,:) = theta.*(2*theta-1);
    N(:,5,:) = 4*xi.*(1-xi-eta-theta);
    N(:,6,:) = 4*xi.*eta;
    N(:,7,:) = 4*eta.*(1-xi-eta-theta);
    N(:,8,:) = 4*theta.*(1-xi-eta-theta);
    N(:,9,:) = 4*xi.*theta;
    N(:,10,:) = 4*eta.*theta;
    if ismydoublend
        N = MYDOUBLEND(N);
    end
else
    % N = inline('[1+xi(1).*(2*xi(1)+4*xi(2)-3)+xi(2).*(2*xi(2)+4*xi(3)-3)+xi(3).*(2*xi(3)+4*xi(1)-3),xi(1).*(2*xi(1)-1),xi(2).*(2*xi(2)-1),xi(3).*(2*xi(3)-1),4*xi(1).*(1-xi(1)-xi(2)-xi(3)),4*xi(1).*xi(2),4*xi(2).*(1-xi(1)-xi(2)-xi(3)),4*xi(3).*(1-xi(1)-xi(2)-xi(3)),4*xi(1).*xi(3),4*xi(2).*xi(3)]','xi');
    N = @(xi) [1+xi(1).*(2*xi(1)+4*xi(2)-3)+xi(2).*(2*xi(2)+4*xi(3)-3)+xi(3).*(2*xi(3)+4*xi(1)-3),...
        xi(1).*(2*xi(1)-1),...
        xi(2).*(2*xi(2)-1),...
        xi(3).*(2*xi(3)-1),...
        4*xi(1).*(1-xi(1)-xi(2)-xi(3)),...
        4*xi(1).*xi(2),...
        4*xi(2).*(1-xi(1)-xi(2)-xi(3)),...
        4*xi(3).*(1-xi(1)-xi(2)-xi(3)),...
        4*xi(1).*xi(3),...
        4*xi(2).*xi(3)];
end
