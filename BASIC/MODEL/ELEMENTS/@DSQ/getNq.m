function N = getNq(elem,xi)
% function N = getNq(elem,xi)

xi = double(xi);
N = zeros([size(xi,1),8,sizeND(xi)]);
eta = xi(:,2,:);
xi = xi(:,1,:); 

N(:,1,:) = -1/4.*(1-xi).*(1-eta).*(1+xi+eta);
N(:,2,:) = -1/4.*(1+xi).*(1-eta).*(1-xi+eta);
N(:,3,:) = 1/4.*(1+xi).*(1+eta).*(-1+xi+eta);
N(:,4,:) = -1/4.*(1-xi).*(1+eta).*(1+xi-eta);
N(:,5,:) = 1/2.*(1-xi.^2).*(1-eta);
N(:,6,:) = 1/2.*(1-eta.^2).*(1+xi);
N(:,7,:) = 1/2.*(1-xi.^2).*(1+eta);
N(:,8,:) = 1/2.*(1-eta.^2).*(1-xi);
