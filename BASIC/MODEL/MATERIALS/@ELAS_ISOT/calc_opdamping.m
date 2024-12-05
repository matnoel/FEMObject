function D = calc_opdamping(mat,elem,xnode,xgauss)
% function D = calc_opdamping(mat,elem,xnode,xgauss)

mu = evalparam(mat,'MU',elem,xnode,xgauss); % damping factor

switch getdim(elem)
    case 1
        S = evalparam(mat,'S',elem,xnode,xgauss); % section
        D = mu*S; % damping operator
    case 2
        e = evalparam(mat,'DIM3',elem,xnode,xgauss); % thickness
        D = e*mu; % damping operator
    case 3
        D = mu; % damping operator
end