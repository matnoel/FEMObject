function D = calc_opmasse(mat,elem,xnode,xgauss)
% function D = calc_opmasse(mat,elem,xnode,xgauss)

rho = evalparam(mat,'RHO',elem,xnode,xgauss); % mass density
S = evalparam(mat,'S',elem,xnode,xgauss); % cross-section area

switch getindim(elem)
    case 2
        D = rho*S*diag([1,1,0]); % mass operator
    case 3
        D = rho*S*diag([1,1,1,0,0,0]); % mass operator
end
