function [se,B] = sigma(mat,elem,xnode,xgauss,qe,varargin)
% function [se,B] = sigma(mat,elem,xnode,xgauss,qe,varargin)

D = calc_opmat(mat,elem,xnode,xgauss);
B = calc_B(elem,xnode,xgauss);

se = D*(B*qe);

if ischarin('local',varargin)
    e = evalparam(mat,'DIM3',elem,xnode,xgauss); % thickness
    Se = se;
    se = zerosND([6,1,sizeND(Se)]);
    % sxx = Nxx/e + Mxx*12/e^3*z with z = -e/2
    % syy = Nyy/e + Myy*12/e^3*z with z = -e/2
    % szz = 0
    % sxy = Nxy/e + Mxy*12/e^3*z with z = -e/2
    se(1) = Se(1)/e - Se(4)*12/e^3*e/2;
    se(2) = Se(2)/e - Se(5)*12/e^3*e/2;
    se(3) = 0;
    se(6) = Se(3)/e - Se(6)*12/e^3*e/2;
    switch class(elem)
        case {'DKT','DKQ'} % Kirchhoff-Love (classical) plate theory
            % syz = 0
            % sxz = 0
            se(4:5) = 0;
        case {'DST','DSQ','COQ4'} % Reissner-Mindlin (first-order shear) plate theory
            % syz = Qy/e
            % sxz = Qx/e
            se(4) = Se(8)/e;
            se(5) = Se(7)/e;
    end
    Pe = calc_Pmat(getsyscoord(elem));
    % [xx yy zz yz xz xy]
    se = Pe'*se;
    
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
