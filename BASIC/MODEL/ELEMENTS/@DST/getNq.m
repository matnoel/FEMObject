function N = getNq(elem,xi)
% function N = getNq(elem,xi)

xi = double(xi);
N = zeros([size(xi,1),6,sizeND(xi)]);
eta = xi(:,2,:);
xi = xi(:,1,:); 

N(:,1,:) = 2*(1-xi-eta).*(1/2-xi-eta);
N(:,2,:) = xi.*(2*xi-1);
N(:,3,:) = eta.*(2*eta-1);
N(:,4,:) = 4*xi.*(1-xi-eta);
N(:,5,:) = 4*xi.*eta;
N(:,6,:) = 4*eta.*(1-xi-eta);
