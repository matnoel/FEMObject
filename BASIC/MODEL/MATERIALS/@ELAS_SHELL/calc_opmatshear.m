function D = calc_opmatshear(mat,elem,xnode,xgauss)
% function D = calc_opmatshear(mat,elem,xnode,xgauss)

e = evalparam(mat,'DIM3',elem,xnode,xgauss); % thickness
E = evalparam(mat,'E',elem,xnode,xgauss); % Young modulus
nu = evalparam(mat,'NU',elem,xnode,xgauss); % Poisson ratio

% isotropic material
G = E/(1+nu)/2; % shear modulus

k = getparam(mat,'k'); % shear correction factor

alpha = k*e*G; % shear stiffness (or shear rigidity)
D = alpha*eye(2); % shear behaviour
