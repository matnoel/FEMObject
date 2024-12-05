function DDNv = getDDNv(elem,xi)
% function DDNv = getDDNv(elem,xi)

if nargin==2
    DDNv = zeros([1,4,sizeND(xi)]);
    DDNv(1,1,:) = 3/2*xi;
    DDNv(1,2,:) = 1/4*(3*xi-1);
    DDNv(1,3,:) = -3/2*xi;
    DDNv(1,4,:) = 1/4*(3*xi+1);
else
    % DDNv = inline('[3/2*xi,1/4*(3*xi-1),-3/2*xi,1/4*(3*xi+1)]','xi');
    DDNv = @(xi) [3/2*xi,1/4*(3*xi-1),-3/2*xi,1/4*(3*xi+1)];
end

DDNv = rescale(DDNv,elem);
