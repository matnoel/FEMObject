function DN = getDNq(elem,xi)
% function DN = getDNq(elem,xi)

xi = double(xi);
DN = zeros([2,6,sizeND(xi)]);
eta = xi(:,2,:);
xi = xi(:,1,:); 

DN(1,1,:) = 4*(xi+eta-3/4);
DN(1,2,:) = 4*xi-1;
DN(1,3,:) = 0*eta;
DN(1,4,:) = 4.*(1-2*xi-eta);
DN(1,5,:) = 4*eta;
DN(1,6,:) = -4*eta;

DN(2,1,:) = 4*(xi+eta-3/4);
DN(2,2,:) = 0*xi;
DN(2,3,:) = 4*eta-1;
DN(2,4,:) = -4*xi;
DN(2,5,:) = 4*xi;
DN(2,6,:) = 4*(1-xi-2*eta);
