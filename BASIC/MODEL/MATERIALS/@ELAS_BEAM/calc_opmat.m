function D = calc_opmat(mat,elem,xnode,xgauss)
% function D = calc_opmat(mat,elem,xnode,xgauss)

E = evalparam(mat,'E',elem,xnode,xgauss); % Young modulus
S = evalparam(mat,'S',elem,xnode,xgauss); % cross-section area

switch getindim(elem)
    case 2
        I = evalparam(mat,'IZ',elem,xnode,xgauss); % planar second moment of area (or planar area moment of inertia)
        D = [E*S,0;...
            0,E*I]; % stiffness operator
        
    case 3
        nu = evalparam(mat,'NU',elem,xnode,xgauss); % Poisson ratio
        G = E/(1+nu)/2; % shear modulus
        IY = evalparam(mat,'IY',elem,xnode,xgauss); % planar second moment of area (or planar area moment of inertia)
        IZ = evalparam(mat,'IZ',elem,xnode,xgauss); % planar second moment of area (or planar area moment of inertia)
        IX = evalparam(mat,'IX',elem,xnode,xgauss); % polar second moment of area (or polar area moment of inertia)
        
        D = [E*S,0,0,0;...
            0,G*IX,0,0;...
            0,0,E*IY,0;...
            0,0,0,E*IZ]; % stiffness operator
end

