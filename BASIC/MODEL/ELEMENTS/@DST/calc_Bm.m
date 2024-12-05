function [B,detJ] = calc_Bm(elem,xnode,xgauss)
% function [B,detJ] = calc_Bm(elem,xnode,xgauss)

nbelem = getnbelem(elem);

[detJ,J,Ji,xnode,DNlocal] = calc_detJ(elem,xnode,xgauss);
DN = Ji*DNlocal;

% MEMBRANE (UX,UY)
B = zerosND([3,3*2,sizeND(DN)]);
B(1,1:2:end) = DN(1,:);
B(2,2:2:end) = DN(2,:);
B(3,2:2:end) = DN(1,:);
B(3,1:2:end) = DN(2,:);
