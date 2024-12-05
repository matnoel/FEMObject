function DDN = getDDNq(elem,xi)
% function DDN = getDDNq(elem,xi)

xi = double(xi);
DDN = zeros([3,6,sizeND(xi)]);

DDN(1,1,:) = 4;
DDN(1,2,:) = 4;
DDN(1,3,:) = 0;
DDN(1,4,:) = -8;
DDN(1,5,:) = 0;
DDN(1,6,:) = 0;

DDN(2,1,:) = 4;
DDN(2,2,:) = 0;
DDN(2,3,:) = 4;
DDN(2,4,:) = 0;
DDN(2,5,:) = 0;
DDN(2,6,:) = -8;

DDN(3,1,:) = 4;
DDN(3,2,:) = 0;
DDN(3,3,:) = 0;
DDN(3,4,:) = -4;
DDN(3,5,:) = 4;
DDN(3,6,:) = -4;
