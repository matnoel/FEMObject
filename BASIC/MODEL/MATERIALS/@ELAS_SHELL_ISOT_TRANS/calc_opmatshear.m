function D = calc_opmatshear(mat,elem,xnode,xgauss)
% function D = calc_opmatshear(mat,elem,xnode,xgauss)

GL = evalparam(mat,'GL',elem,xnode,xgauss); % longitudinal shear modulus
k = getparam(mat,'k'); % shear correction factor
e = evalparam(mat,'DIM3',elem,xnode,xgauss); % thickness

alpha = k*e*GL; % shear stiffness (or shear rigidity)
D = alpha*eye(2); % shear behaviour
