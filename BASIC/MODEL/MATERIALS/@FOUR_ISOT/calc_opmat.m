function k = calc_opmat(mat,elem,xnode,xgauss)
% function k = calc_opmat(mat,elem,xnode,xgauss)

if nargin<=2
    xnode = [];
    xgauss = [];
end

k = evalparam(mat,'k',elem,xnode,xgauss);

switch getdim(elem)
    case 1
        if isparam(mat,'S')
            S = evalparam(mat,'S',elem,xnode,xgauss); % cross-section area
            k = k*S;
        end
    case 2
        if isparam(mat,'DIM3')
            e = evalparam(mat,'DIM3',elem,xnode,xgauss); % thickness
            k = k*e;
        end
end
