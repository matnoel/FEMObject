function [B,detJ] = calc_Bd(elem,xnode,xgauss)
% function [B,detJ] = calc_Bd(elem,xnode,xgauss)

nbelem = getnbelem(elem);

N = MYDOUBLEND(getN(elem,xgauss));

% DRILLING RZ
B = zerosND([1,4,sizeND(N)]);
B(1,:) = N(1,:);
