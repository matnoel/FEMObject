function N = getN(elem,xi)
% function N = getN(elem,xi)

if nargin==2
    s = sizeND(xi);
    xi = double(xi);
    N = zeros([size(xi,1),8,s]);
    
    fact1 = repmat([-1,1,1,-1,-1,1,1,-1],[size(xi,1),1,prod(s)]);
    fact2 = repmat([-1,-1,1,1,-1,-1,1,1],[size(xi,1),1,prod(s)]);
    fact3 = repmat([-1,-1,-1,-1,1,1,1,1],[size(xi,1),1,prod(s)]);
    xi1 = repmat(xi(:,1,:),[1,8,1]);
    xi2 = repmat(xi(:,2,:),[1,8,1]);
    xi3 = repmat(xi(:,3,:),[1,8,1]);
    
    N(:,:,:) = (1+xi1.*fact1).*(1+xi2.*fact2).*(1+xi3.*fact3)/8;
else
    % N = inline('1/8*[(1+xi(:,1)*[-1,1,1,-1,-1,1,1,-1]).*(1+xi(:,2)*[-1,-1,1,1,-1,-1,1,1]).*(1+xi(:,3)*[-1,-1,-1,-1,1,1,1,1])]','xi');
    N = @(xi) 1/8*[(1+xi(:,1)*[-1,1,1,-1,-1,1,1,-1]).*(1+xi(:,2)*[-1,-1,1,1,-1,-1,1,1]).*(1+xi(:,3)*[-1,-1,-1,-1,1,1,1,1])];
end
