function ke = rigitangexpect(mat,elem,xnode,xgauss,qe,l1,l2,varargin)
% function ke = rigitangexpect(mat,elem,xnode,xgauss,qe,l1,l2,varargin)

N = calc_N(elem,xnode,xgauss);
u = N*qe;
k = evalparampc(mat,'k',PC,elem,xnode,xgauss);
B = calc_B(elem,xnode,xgauss);
k = expect(k,l1,l2);

switch getparam(mat,'formulation')
    case 1
        if ~israndom(u)
            ke = B'*k*B - B'*(1/2-u)*N;
        else
            ke = B'*k*B - B'*(1/2-expect(u,l1,l2))*N;
        end
        
    case 2
        gradu = B*qe;
        alpha = getparam(mat,'alpha');
        Ntilde = N + alpha*B;
        if ~israndom(u)
            ke = B'*k*B + Ntilde'*(1/2*B - gradu*N - u*B);
        else
            ke = B'*k*B + Ntilde'*(1/2*B - expect(gradu,l1,l2)*N - expect(u,l1,l2)*B);
        end
        
end
