function [B,detJ] = calc_Bb(elem,xnode,xgauss)
% function [B,detJ] = calc_Bb(elem,xnode,xgauss)

nbelem = getnbelem(elem);

[detJ,J,Ji,xnode,DNlocal] = calc_detJ(elem,xnode,xgauss);
DN = Ji*DNlocal;

% BENDING (RX,RY)
B = zerosND([3,4*2,sizeND(DN)]);
B(1,2:2:end) = DN(1,:);
B(2,1:2:end) = -DN(2,:);
B(3,1:2:end) = -DN(1,:);
B(3,2:2:end) = DN(2,:);
