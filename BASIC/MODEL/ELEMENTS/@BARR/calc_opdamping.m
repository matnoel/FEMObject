function D = calc_opdamping(elem,mat,xnode,xgauss)
% function D = calc_opdamping(elem,mat,xnode,xgauss)

switch class(mat)
    case {'ELAS_ISOT','ELAS_ISOT_TRANS'}
        mu = getparam(mat,'MU') ;
        S = getparam(mat,'S') ;
        D = mu*S*speye(getindim(elem));
end