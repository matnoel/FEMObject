function ke = calc_rigielem(elem,node,mat,gauss,varargin)
% function ke = calc_rigielem(elem,node,mat,gauss)
ke=sparse(elem.nbddl,elem.nbddl);

node=getcoord(node);
for l=1:gauss.nbgauss
   [B,detJ]=calc_B(elem,node,gauss.coord(l,:));
   D=calc_opmat(elem,mat,gauss.coord(l,:));  
   ke=ke+gauss.w(l)*abs(detJ)*B'*D*B;  
end
