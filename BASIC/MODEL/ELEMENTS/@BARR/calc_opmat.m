function D = calc_opmat(elem,mat,xnode,xgauss)
% function D = calc_opmat(elem,mat,xnode,xgauss)

switch class(mat)
    case {'ELAS_ISOT'}
        E = evalparam(elem,mat,'E',xnode,xgauss); % Young modulus
        S = evalparam(elem,mat,'S',xnode,xgauss); % cross-section area
        D = E*S; % stiffness operator
end