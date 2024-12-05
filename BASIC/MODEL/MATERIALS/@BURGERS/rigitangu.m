function keU = rigitangu(mat,elem,xnode,xgauss,qe,Ue,varargin)
% function keU = rigitangu(mat,elem,xnode,xgauss,qe,Ue,varargin)

N = calc_N(elem,xnode,xgauss);
u = N*qe;
k = evalparampc(mat,'k',[],elem,xnode,xgauss);
B = calc_B(elem,xnode,xgauss);
BU = B*Ue;

switch getparam(mat,'formulation')
    case 1
        keU = B'*k*BU - B'*(1/2-u)*NU;
        
    case 2
        gradu = B*qe;
        alpha = getparam(mat,'alpha');
        Ntilde = N + alpha*B;
        keU = B'*k*BU + Ntilde'*(1/2*BU - gradu*NU - u*BU);
        
end
