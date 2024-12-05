function D = calc_opmat(mat,elem,xnode,xgauss)
% function D = calc_opmat(mat,elem,xnode,xgauss)

E = evalparam(mat,'E',elem,xnode,xgauss); % Young modulus
nu = evalparam(mat,'NU',elem,xnode,xgauss); % Poisson ratio
e = evalparam(mat,'DIM3',elem,xnode,xgauss); % thickness

% isotropic material
G = E/(1+nu)/2; % shear modulus

k = getparam(mat,'k'); % shear correction factor

a = e*E/(1-nu^2); % extensional stiffness (or membrane rigidity)
d = e^3*E/(1-nu^2)/12; % bending stiffness (or flexural rigidity)
alpha = k*e*G; % shear stiffness (or shear rigidity)

C1 = a*[1,nu,0;nu,1,0;0,0,(1-nu)/2]; % membrane behaviour
C2 = d*[1,nu,0;nu,1,0;0,0,(1-nu)/2]; % bending behaviour
C3 = alpha*eye(2); % shear behaviour
drill = 1e-10*E*e.^1; % drilling behaviour

switch class(elem)
   case {'DKT','DKQ'} % Kirchhoff-Love (classical) plate theory
       D = [C1,zeros(3,3),zeros(3,1);...
           zeros(3,3),C2,zeros(3,1);...
           zeros(1,3),zeros(1,3),drill]; % stiffness operator
    case {'DST','DSQ','COQ4'} % Reissner-Mindlin (first-order shear) plate theory
        D = [C1,zeros(3,3),zeros(3,2),zeros(3,1);...
            zeros(3,3),C2,zeros(3,2),zeros(3,1);...
            zeros(2,3),zeros(2,3),C3,zeros(2,1);...
            zeros(1,3),zeros(1,3),zeros(1,2),drill]; % stiffness operator
    otherwise
        error(['Wrong elem type ' class(elem) ' for material type ' class(mat)])
        
end
