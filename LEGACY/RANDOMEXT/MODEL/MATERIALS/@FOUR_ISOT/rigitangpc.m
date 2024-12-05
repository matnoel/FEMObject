function ke = rigitangpc(mat,elem,xnode,xgauss,qe,PC,varargin)
% function ke = rigitangpc(mat,elem,xnode,xgauss,qe,PC,varargin)

k = evalparampc(mat,'k',PC,elem,xnode,xgauss);
B = calc_B(elem,xnode,xgauss);
ke = B'*k*B;

if mat.b || mat.k2 || mat.r
    N = calc_N(elem,xnode,xgauss);
    
    if mat.b
        b = evalparampc(mat,'b',PC,elem,xnode,xgauss);
        keb = (B'*b)*N;
        ke = ke + keb;
    end
    
    if mat.k2
        k2 = evalparampc(mat,'k2',PC,elem,xnode,xgauss);
        % k2 = expand(k2);
        u = N*qe;
        du = B*qe;
        
        kek2 = B'*((k2.*u.*u)*B+(2.*k2.*du)*N);
        
        ke = ke + kek2;
    end
    
    if mat.r
        r = evalparampc(mat,'r',PC,elem,xnode,xgauss);
        % r = expand(r);
        u = N*qe;
        
        ker = N'*((3.*r.*u.*u)*N);
        
        ke = ke + ker;
    end
end
