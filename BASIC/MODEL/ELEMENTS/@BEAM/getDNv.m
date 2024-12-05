function DNv = getDNv(elem,xi)
% function DNv = getDNv(elem,xi)

if nargin==2
    DNv = zeros([1,4,sizeND(xi)]);
    DNv(1,1,:) = -3/4*(1-xi).*(1+xi);
    DNv(1,2,:) = -1/8*(1-xi).*(1+3*xi);
    DNv(1,3,:) = 3/4*(1-xi).*(1+xi);
    DNv(1,4,:) = 1/8*(1+xi).*(3*xi-1);
else
    % DNv = inline('[-3/4*(1-xi)*(1+xi),-1/8*(1-xi)*(1+3*xi),3/4*(1-xi)*(1+xi),1/8*(1+xi)*(3*xi-1)]','xi');
    DNv = @(xi) [-3/4*(1-xi)*(1+xi),-1/8*(1-xi)*(1+3*xi),3/4*(1-xi)*(1+xi),1/8*(1+xi)*(3*xi-1)];
end

DNv = rescale(DNv,elem);
