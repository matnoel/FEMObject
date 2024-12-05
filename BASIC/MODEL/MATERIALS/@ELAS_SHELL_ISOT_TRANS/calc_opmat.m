function D = calc_opmat(mat,elem,xnode,xgauss)
% function D = calc_opmat(mat,elem,xnode,xgauss)

ET = evalparam(mat,'ET',elem,xnode,xgauss); % transverse Young modulus
nuT = evalparam(mat,'NUT',elem,xnode,xgauss); % transverse Poisson ratio
GL = evalparam(mat,'GL',elem,xnode,xgauss); % longitudinal shear modulus
e = evalparam(mat,'DIM3',elem,xnode,xgauss); % thickness
k = getparam(mat,'k'); % shear correction factor

a = e*ET/(1-nuT^2); % extensional stiffness (or membrane rigidity)
d = e^3*ET/(1-nuT^2)/12; % bending stiffness (or flexural rigidity)
alpha = k*e*GL; % shear stiffness (or shear rigidity)

C1 = a*[1,nuT,0;nuT,1,0;0,0,(1-nuT)/2]; % membrane behaviour
C2 = d*[1,nuT,0;nuT,1,0;0,0,(1-nuT)/2]; % bending behaviour
C3 = alpha*eye(2); % shear behaviour
drill = 1e-10*ET*e.^1; % drilling behaviour

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
