function [B,detJ] = calc_B(elem,xnode,xgauss)
% function [B,detJ] = calc_B(elem,xnode,xgauss)

[DN,detJ] = calc_DN(elem,xnode,xgauss);
P = calc_Plocal(elem,xnode);
B = DN*P;
