function [x,Nuni] = calc_x(elem,xnode,xi)

xnode = calc_xnode(elem,xnode);

Nuni = getN(elem,xi);

x=Nuni*xnode;
