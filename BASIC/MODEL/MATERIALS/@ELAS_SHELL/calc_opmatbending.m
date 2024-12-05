function D = calc_opmatbending(mat,elem,xnode,xgauss)
% function D = calc_opmatbending(mat,elem,xnode,xgauss)

E = evalparam(mat,'E',elem,xnode,xgauss); % Young modulus
nu = evalparam(mat,'NU',elem,xnode,xgauss); % Poisson ratio
e = evalparam(mat,'DIM3',elem,xnode,xgauss); % thickness

d = e^3*E/(1-nu^2)/12; % bending stiffness (or flexural rigidity)
D = d*[1,nu,0;nu,1,0;0,0,(1-nu)/2]; % bending behaviour
