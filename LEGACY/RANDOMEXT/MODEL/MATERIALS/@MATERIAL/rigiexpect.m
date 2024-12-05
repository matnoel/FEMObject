function ke = rigiexpect(mat,elem,xnode,xi,a,b,varargin)
% function ke = rigiexpect(mat,elem,xnode,xi,a,b,varargin)

D = calc_opmatpc(mat,elem,xnode,xi,[]);
D = expect(D,a,b);
B = calc_B(elem,xnode,xi);
ke = B'*D*B;

