function D = calc_opmasse(mat,elem,xnode,xgauss)
% function D = calc_opmasse(mat,elem,xnode,xgauss)

rho = evalparam(mat,'RHO',elem,xnode,xgauss); % mass density
e = evalparam(mat,'DIM3',elem,xnode,xgauss); % thickness

D = rho*e*[eye(3),zeros(3,3);zeros(3,3),zeros(3,3)]; % mass operator
