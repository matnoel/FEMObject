function D = calc_opdamping(mat,elem,xnode,xgauss)
% function D = calc_opdamping(mat,elem,xnode,xgauss)

mu = evalparam(mat,'MU',elem,xnode,xgauss); % damping factor
S = evalparam(mat,'S',elem,xnode,xgauss); % cross-section area

switch getindim(elem)
    case 2
        D = mu*S*diag([1,1,0]); % damping operator
    case 3
        D = mu*S*diag([1,1,1,0,0,0]); % damping operator
end
