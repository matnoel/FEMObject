function ke = rigipc(mat,elem,xnode,xi,PC,varargin)
% function ke = rigipc(mat,elem,xnode,xi,PC,varargin)

D = calc_opmatpc(mat,elem,xnode,xi,PC);
B = calc_B(elem,xnode,xi);
ke = B'*D*B;
