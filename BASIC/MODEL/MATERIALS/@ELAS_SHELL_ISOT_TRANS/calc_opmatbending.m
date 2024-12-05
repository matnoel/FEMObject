function D = calc_opmatbending(mat,elem,xnode,xgauss)
% function D = calc_opmatbending(mat,elem,xnode,xgauss)

ET = evalparam(mat,'ET',elem,xnode,xgauss); % Young modulus
e = evalparam(mat,'DIM3',elem,xnode,xgauss); % thickness
nuT = evalparam(mat,'NUT',elem,xnode,xgauss); % Poisson ratio

d = e^3*ET/(1-nuT^2)/12; % bending stiffness (or flexural rigidity)
D = d*[1,nuT,0;nuT,1,0;0,0,(1-nuT)/2]; % bending behaviour
