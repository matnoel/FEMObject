function ke = rigitang(mat,elem,xnode,xgauss,qe,varargin)
% function ke = rigitang(mat,elem,xnode,xgauss,qe,varargin)

N = calc_N(elem,xnode,xgauss);
u = N*qe;
k = evalparam(mat,'k',elem,xnode,xgauss);
B = calc_B(elem,xnode,xgauss);

switch getparam(mat,'formulation')
    case 1
        ke = B'*k*B - B'*(1/2-u)*N;
        
    case 2
        gradu = B*qe;
        alpha = getparam(mat,'alpha');
        Ntilde = N + alpha*B;
        ke = B'*k*B + Ntilde'*(1/2*B - gradu*N - u*B);
        
end
