function D = calc_opmassepc(mat,elem,xnode,xgauss,PC)
% function D = calc_opmassepc(mat,elem,xnode,xgauss,PC)

rho = evalparampc(mat,'RHO',PC,elem,xnode,xgauss);
switch getdim(elem)
    case 1
        S = evalparampc(mat,'S',PC,elem,xnode,xgauss);
        D = rho*S;
    case 2
        e = evalparampc(mat,'DIM3',PC,elem,xnode,xgauss);
        D = e*rho;
    case 3
        D = rho;
end
