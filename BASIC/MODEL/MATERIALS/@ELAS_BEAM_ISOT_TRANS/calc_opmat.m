function D = calc_opmat(mat,elem,xnode,xgauss)
% function D = calc_opmat(mat,elem,xnode,xgauss)

EL = evalparam(mat,'EL',elem,xnode,xgauss); % longitudinal Young modulus
S = evalparam(mat,'S',elem,xnode,xgauss); % cross-section area

switch getindim(elem)
    case 2
        I = evalparam(mat,'IZ',elem,xnode,xgauss); % planar second moment of area (or planar area moment of inertia)
        D = [EL*S,0;...
            0,EL*I]; % stiffness operator
        
    case 3
        nuT = evalparam(mat,'NUT',elem,xnode,xgauss); % transverse Poisson ratio
        GT = EL/(1+nuT)/2; % transverse shear modulus
        GL = evalparam(mat,'GL',elem,xnode,xgauss); % longitudinal shear modulus
        IY = evalparam(mat,'IY',elem,xnode,xgauss); % planar second moment of area (or planar area moment of inertia)
        IZ = evalparam(mat,'IZ',elem,xnode,xgauss); % planar second moment of area (or planar area moment of inertia)
        IX = evalparam(mat,'IX',elem,xnode,xgauss); % polar second moment of area (or polar area moment of inertia)
        
        D = [EL*S,0,0,0;...
            0,GL*IX,0,0;...
            0,0,EL*IY,0;...
            0,0,0,EL*IZ]; % stiffness operator
end

