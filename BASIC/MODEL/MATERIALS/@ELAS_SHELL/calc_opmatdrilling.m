function D = calc_opmatdrilling(mat,elem,xnode,xgauss)
% function D = calc_opmatdrilling(mat,elem,xnode,xgauss)

E = evalparam(mat,'E',elem,xnode,xgauss); % Young modulus
e = evalparam(mat,'DIM3',elem,xnode,xgauss); % thickness

D = 1e-10*E*e^1; % drilling behaviour
