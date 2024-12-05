function P = calc_Plocal(elem,xnode)
% function P = calc_Plocal(elem,xnode)

R = calc_Rlocal(elem);
z = zerosND(size(R));
P = [R,z;z,R];
