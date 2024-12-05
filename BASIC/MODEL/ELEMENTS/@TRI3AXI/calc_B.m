function [B,detJ,x] = calc_B(elem,xnode,xgauss)
% function [B,detJ,x] = calc_B(elem,xnode,xgauss)

if nargin<3 || isempty(xgauss)
    xgauss = [1/3,1/3];
end
BT = calc_B(elem.TRI3,xnode,xgauss);
N = calc_N(elem.TRI3,xnode,xgauss);
x = calc_x(elem,xnode,xgauss);
B = zerosND([4,6,sizeND(N)]);
B([1,2,4],:) = BT;
B(3,1:2:end) = N(1,1:2:end)./x(1,1);

if nargin>1
    detJ = calc_detJ(elem,xnode,xgauss);
end
