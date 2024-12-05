function D = calc_opmatdrilling(mat,elem,xnode,xgauss)
% function D = calc_opmatdrilling(mat,elem,xnode,xgauss)

ET = evalparam(mat,'ET',elem,xnode,xgauss); % transverse Young modulus
e = evalparam(mat,'DIM3',elem,xnode,xgauss); % thickness

D = 1e-10*ET*e^1; % drilling behaviour
