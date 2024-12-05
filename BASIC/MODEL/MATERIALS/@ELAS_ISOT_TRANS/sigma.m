function [se,B] = sigma(mat,elem,xnode,xgauss,qe,varargin)
% function [se,B] = sigma(mat,elem,xnode,xgauss,qe,varargin)

B = calc_B(elem,xnode,xgauss);
D = calc_opmat(mat,elem,xnode,xgauss);

se = D*(B*qe);

if ischarin('local',varargin)
    switch getdim(elem)
        case 1
            S = evalparam(mat,'S',elem,xnode,xgauss); % cross-section area
            se = se/S;
        case 2
            e = evalparam(mat,'DIM3',elem,xnode,xgauss); % thickness
            se = se/e;
    end
end
