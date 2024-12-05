function DDN = getDDNq(elem,xi)
% function DDN = getDDNq(elem,xi)

xi = double(xi);
DDN = zeros([3,8,sizeND(xi)]);
eta = xi(:,2,:);
xi = xi(:,1,:); 

DDN(1,1,:) = 1/2.*(1-eta);
DDN(1,2,:) = 1/2.*(1-eta);
DDN(1,3,:) = 1/2.*(1+eta);
DDN(1,4,:) = 1/2.*(1+eta);
DDN(1,5,:) = eta-1;
DDN(1,6,:) = 0;
DDN(1,7,:) = -(1+eta);
DDN(1,8,:) = 0;

DDN(2,1,:) = 1/2.*(1-xi);
DDN(2,2,:) = 1/2.*(1+xi);
DDN(2,3,:) = 1/2.*(1+xi);
DDN(2,4,:) = 1/2.*(1-xi);
DDN(2,5,:) = 0;
DDN(2,6,:) = -(1+xi);
DDN(2,7,:) = 0;
DDN(2,8,:) = xi-1;

DDN(3,1,:) = 1/4.*(1-2.*eta-2.*xi);
DDN(3,2,:) = 1/4.*(2.*eta-2.*xi-1);
DDN(3,3,:) = 1/4.*(1+2.*eta+2.*xi);
DDN(3,4,:) = 1/4.*(2.*xi-2.*eta-1);
DDN(3,5,:) = xi;
DDN(3,6,:) = -eta;
DDN(3,7,:) = -xi;
DDN(3,8,:) = eta;
