function [se,B] = sigma(mat,elem,xnode,xgauss,qe,varargin)
% function [se,B] = sigma(mat,elem,xnode,xgauss,qe,varargin)

D = calc_opmat(mat,elem,xnode,xgauss);
B = calc_B(elem,xnode,xgauss);

se = D*(B*qe);

if ischarin('local',varargin)
    Se = se;
    S = evalparam(mat,'S',elem,xnode,xgauss); % cross-section area
    switch getindim(elem)
        case 2
            se = zerosND([3,1,sizeND(Se)]);
            I = evalparam(mat,'IZ',elem,xnode,xgauss); % planar second moment of area (or planar area moment of inertia)
            y = getcharin('y',varargin,sqrt(S));
            % sxx = N/S - Mz/Iz*y
            % syy = 0
            % sxy = Ty/S
            se(1) = Se(1)/S - Se(2)/I*y;
            se(2) = 0;
            se(3) = 0;
        case 3
            se = zerosND([6,1,sizeND(Se)]);
            IY = evalparam(mat,'IY',elem,xnode,xgauss); % planar second moment of area (or planar area moment of inertia)
            IZ = evalparam(mat,'IZ',elem,xnode,xgauss); % planar second moment of area (or planar area moment of inertia)
            IX = evalparam(mat,'IX',elem,xnode,xgauss); % polar second moment of area (or polar area moment of inertia)
            y = getcharin('y',varargin,sqrt(S));
            z = getcharin('z',varargin,sqrt(S));
            % sxx = N/S + My/Iy*z - Mz/Iz*y
            % syy = 0
            % szz = 0
            % syz = 0
            % sxz = Tz/S + Mx/Ix*y
            % sxy = Ty/S - Mx/Ix*z
            se(1) = Se(1)/S + Se(3)/IY*z - Se(4)/IZ*y;
            se(2) = 0;
            se(3) = 0;
            se(4) = 0;
            se(5) = Se(2)/IX*y;
            se(6) = -Se(2)/IX*z;
    end
    Pe = calc_Pmat(getsyscoord(elem));
    % [xx yy xy] in 2D
    % [xx yy zz yz xz xy] in 3D
    se = Pe'*se;
    
    if getindim(elem)==3
        P = [1,0,0,0,0,0;...
            0,1,0,0,0,0;...
            0,0,1,0,0,0;...
            0,0,0,0,0,1;...
            0,0,0,0,1,0;...
            0,0,0,1,0,0];
        % [xx yy zz xy xz yz]
        se = P'*se;
    end
end
 
end
