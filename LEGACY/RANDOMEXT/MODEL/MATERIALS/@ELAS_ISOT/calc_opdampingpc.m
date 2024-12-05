function D = calc_opdampingpc(mat,elem,xnode,xgauss,PC)
% function D = calc_opdampingpc(mat,elem,xnode,xgauss,PC)

mu = evalparampc(mat,'MU',PC,elem,xnode,xgauss);
switch getdim(elem)
    case 1
        S = evalparampc(mat,'S',PC,elem,xnode,xgauss);
        D = mu*S;
    case 2
        e = evalparampc(mat,'DIM3',PC,elem,xnode,xgauss);
        D = e*mu;
    case 3
        D = mu;
end
