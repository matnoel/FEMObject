function [B,detJ] = calc_Bs(elem,xnode,xgauss)
% function [B,detJ] = calc_Bs(elem,xnode,xgauss)

nbelem = getnbelem(elem);

N = MYDOUBLEND(getN(elem,xgauss));
[detJ,J,Ji,xnode,DNlocal] = calc_detJ(elem,xnode,xgauss);
DN = Ji*DNlocal;

% SHEAR (UZ,RX,RY)
B = zerosND([2,4*3,sizeND(DN)]);
B(1,1:3:end) = DN(1,:); % UZ'+RY
B(1,3:3:end) = N(1,:);
B(2,1:3:end) = DN(2,:); % UZ'-RX
B(2,2:3:end) = -N(1,:);
