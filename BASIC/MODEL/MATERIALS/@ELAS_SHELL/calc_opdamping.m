function D = calc_opdamping(mat,elem,xnode,xgauss)
% function D = calc_opdamping(mat,elem,xnode,xgauss)

mu = evalparam(mat,'MU',elem,xnode,xgauss); % damping factor
e = evalparam(mat,'DIM3',elem,xnode,xgauss); % thickness

D = e*mu*[eye(3),zeros(3,3);zeros(3,3),zeros(3,3)]; % damping operator
