function [se,B] = sigmapc(mat,elem,xnode,xgauss,qe,PC,varargin)
% function [se,B] = sigmapc(mat,elem,xnode,xgauss,qe,PC,varargin)

k = evalparampc(mat,'k',PC,elem,xnode,xgauss);
B = calc_B(elem,xnode,xgauss);
gradu = B*qe;

se = k.*gradu;

if mat.b || mat.k2 || mat.r
    N = calc_N(elem,xnode,xgauss);
    
    warning('mal calcule')
    if mat.b
        b = evalparampc(mat,'b',PC,elem,xnode,xgauss);
        se = se + b'*gradu;
    end
    
    if mat.k2
        u = N*qe;
        k2 = evalparampc(mat,'k2',PC,elem,xnode,xgauss);
        se = se +(u.*u).*(k2.*gradu);
    end
    
    if mat.r
        u = N*qe;
        r = evalparampc(mat,'r',PC,elem,xnode,xgauss);
        se = se +(u.*u).*(r.*u);
    end
end


