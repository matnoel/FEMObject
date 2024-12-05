function fe = fintpc(mat,elem,xnode,xgauss,qe,PC,varargin)
% function fe = fintpc(mat,elem,xnode,xgauss,qe,PC,varargin)

k = evalparampc(mat,'k',PC,elem,xnode,xgauss);
B = calc_B(elem,xnode,xgauss);
gradu = B*qe;

fe = B'*(k.*gradu);

if mat.b || mat.k2 || mat.r
    N = calc_N(elem,xnode,xgauss);
    
    if mat.b
        b = evalparampc(mat,'b',PC,elem,xnode,xgauss);
        fe = fe + N'*(b'*gradu);
    end
    
    if mat.k2
        u = N*qe;
        k2 = evalparampc(mat,'k2',PC,elem,xnode,xgauss);
        fe = fe + B'*((u.*u).*(k2.*gradu));
    end
    
    if mat.r
        r = evalparampc(mat,'r',PC,elem,xnode,xgauss);
        fe = fe + N'*((u.*u.*u).*(r.*u)./4);
    end
end


