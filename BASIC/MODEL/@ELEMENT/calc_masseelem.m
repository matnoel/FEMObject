function me = calc_masseelem(elem,node,mat,gauss,varargin)
% function me = calc_masseelem(elem,node,mat,gauss)

me=sparse(elem.nbddl,elem.nbddl);
node=getcoord(node);
for l=1:gauss.nbgauss
   [N,detJ]=calc_N(elem,node,gauss.coord(l,:),varargin{:});
   D=calc_opmasse(elem,mat,gauss.coord(l,:));
   me=me+gauss.w(l)*abs(detJ)*N'*D*N;  
end
