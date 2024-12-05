function D = calc_opmasse(elem,mat,xnode,xgauss)
% function D = calc_opmasse(elem,mat,xnode,xgauss)

switch class(mat)
    case {'ELAS_ISOT'}
        rho = getparam(mat,'RHO') ;
        S = getparam(mat,'S') ;
        D = rho*S*speye(getindim(elem));
end