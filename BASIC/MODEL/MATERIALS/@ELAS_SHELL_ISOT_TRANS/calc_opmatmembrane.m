function D = calc_opmatmembrane(mat,elem,xnode,xgauss)
% function D = calc_opmatmembrane(mat,elem,xnode,xgauss)

ET = evalparam(mat,'ET',elem,xnode,xgauss); % transverse Young modulus
nuT = evalparam(mat,'NUT',elem,xnode,xgauss); % transverse Poisson ratio
e = evalparam(mat,'DIM3',elem,xnode,xgauss); % thickness

a = e*ET/(1-nuT^2); % extensional stiffness (or membrane rigidity)
D = a*[1,nuT,0;nuT,1,0;0,0,(1-nuT)/2]; % membrane behaviour
