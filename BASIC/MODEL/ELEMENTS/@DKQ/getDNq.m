function DN = getDNq(elem,xi)
% function DN = getDNq(elem,xi)

xi = double(xi);
DN = zeros([2,8,sizeND(xi)]);
eta = xi(:,2,:);
xi = xi(:,1,:); 

DN(1,1,:) = 1/4.*(1-eta).*(2.*xi+eta);
DN(1,2,:) = 1/4.*(1-eta).*(2.*xi-eta);
DN(1,3,:) = 1/4.*(1+eta).*(2.*xi+eta);
DN(1,4,:) = 1/4.*(1+eta).*(2.*xi-eta);
DN(1,5,:) = -xi.*(1-eta);
DN(1,6,:) = 1/2.*(1-eta.^2);
DN(1,7,:) = -xi.*(1+eta);
DN(1,8,:) = -1/2.*(1-eta.^2);

DN(2,1,:) = 1/4.*(1-xi).*(2.*eta+xi);
DN(2,2,:) = 1/4.*(1+xi).*(2.*eta-xi);
DN(2,3,:) = 1/4.*(1+xi).*(2.*eta+xi);
DN(2,4,:) = 1/4.*(1-xi).*(2.*eta-xi);
DN(2,5,:) = -1/2.*(1-xi.^2);
DN(2,6,:) = -eta.*(1+xi);
DN(2,7,:) = 1/2.*(1-xi.^2);
DN(2,8,:) = -eta.*(1-xi);
