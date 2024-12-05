function fe = fintpc(mat,elem,xnode,xi,q,PC,varargin)
% function fe = fintpc(mat,elem,xnode,xi,q,PC,varargin)

D = calc_opmatpc(mat,elem,xnode,xi,PC);
B = calc_B(elem,xnode,xi);
error('pas programme')
fe = B'*D*(B*q);

