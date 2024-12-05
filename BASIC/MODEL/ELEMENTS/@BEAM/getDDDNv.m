function DDDNv = getDDDNv(elem,xi)
% function DDDNv = getDDDNv(elem,xi)

if nargin==2
    DDDNv = zeros([1,4,sizeND(xi)]);
    DDDNv(1,1,:) = 3/2;
    DDDNv(1,2,:) = 3/4;
    DDDNv(1,3,:) = -3/2;
    DDDNv(1,4,:) = 3/4;
else
    % DDDNv = inline('[3/2*xi,1/4*(3*xi-1),-3/2*xi,1/4*(3*xi+1)]','xi');
    DDDNv = @(xi) [3/2*xi,1/4*(3*xi-1),-3/2*xi,1/4*(3*xi+1)];
end

DDDNv = rescale(DDDNv,elem);
