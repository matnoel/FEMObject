function c = calc_opmasse(mat,elem,xnode,xgauss)
% function c = calc_opmasse(mat,elem,xnode,xgauss)

c = evalparam(mat,'c',elem,xnode,xgauss);

switch getdim(elem)
    case 1
        if isparam(mat,'S')
            S = evalparam(mat,'S',elem,xnode,xgauss); % cross-section area
            c = c*S;
        end
    case 2
        if isparam(mat,'DIM3')
            e = evalparam(mat,'DIM3',elem,xnode,xgauss); % thickness
            c = c*e;
        end
end