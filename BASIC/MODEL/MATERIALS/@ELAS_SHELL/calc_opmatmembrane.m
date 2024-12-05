function D = calc_opmatmembrane(mat,elem,xnode,xgauss)
% function D = calc_opmatmembrane(mat,elem,xnode,xgauss)

E = evalparam(mat,'E',elem,xnode,xgauss); % Young modulus
e = evalparam(mat,'DIM3',elem,xnode,xgauss); % thickness
nu = evalparam(mat,'NU',elem,xnode,xgauss); % Poisson ratio

a = e*E/(1-nu^2); % extensional stiffness (or membrane rigidity)
D = a*[1,nu,0;nu,1,0;0,0,(1-nu)/2]; % membrane behaviour
