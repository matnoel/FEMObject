function DN = getDN(elem,xi)
% function DN = getDN(elem,xi)

if nargin==2
    xi = double(xi);
    DN = zeros([3,10,sizeND(xi)]);
    theta = xi(:,3,:);
    eta = xi(:,2,:);
    xi = xi(:,1,:);
    
    DN(1,1,:) = 4*(xi+eta+theta-3/4);
    DN(1,2,:) = 4*xi-1;
    DN(1,3,:) = 0*eta;
    DN(1,4,:) = 0*theta;
    DN(1,5,:) = 4.*(1-2*xi-eta-theta);
    DN(1,6,:) = 4*eta;
    DN(1,7,:) = -4*eta;
    DN(1,8,:) = -4*theta;
    DN(1,9,:) = 4*theta;
    DN(1,10,:) = 0*xi;
    
    DN(2,1,:) = 4*(xi+eta+theta-3/4);
    DN(2,2,:) = 0*xi;
    DN(2,3,:) = 4*eta-1;
    DN(2,4,:) = 0*theta;
    DN(2,5,:) = -4*xi;
    DN(2,6,:) = 4*xi;
    DN(2,7,:) = 4.*(1-2*eta-xi-theta);
    DN(2,8,:) = -4*theta;
    DN(2,9,:) = 0*eta;
    DN(2,10,:) = 4*theta;
    
    DN(3,1,:) = 4*(xi+eta+theta-3/4);
    DN(3,2,:) = 0*xi;
    DN(3,3,:) = 0*eta;
    DN(3,4,:) = 4*theta-1;
    DN(3,5,:) = -4*xi;
    DN(3,6,:) = 0*theta;
    DN(3,7,:) = -4*eta;
    DN(3,8,:) = 4.*(1-2*theta-xi-eta);
    DN(3,9,:) = 4*xi;
    DN(3,10,:) = 4*eta;
end
