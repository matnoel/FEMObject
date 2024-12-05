function ke = rigitangexpect(mat,elem,xnode,xgauss,qe,l1,l2,varargin)
% function ke = rigitangexpect(mat,elem,xnode,xgauss,qe,l1,l2,varargin)

k = evalparampc(mat,'k',PC,elem,xnode,xgauss);
B = calc_B(elem,xnode,xgauss);
k = expect(k,l1,l2);
ke = B'*k*B;

if mat.b || mat.k2 || mat.r
    N = calc_N(elem,xnode,xgauss);
    
    if mat.b
        b = evalparampc(mat,'b',PC,elem,xnode,xgauss);
        b = expect(b,l1,l2);
        keb = (B'*b)*N;
        ke = ke + keb;
    end
    
    if mat.k2
        k2 = evalparampc(mat,'k2',PC,elem,xnode,xgauss);
        u = N*qe;
        du = B*qe;
        if ~israndom(u)
            k2 = expect(k2,l1,l2);
            kek2 = B'*((k2.*u.*u)*B+(2.*k2.*du)*N);
        else
            kek2 = B'*(expect(k2.*u.*u,l1,l2)*B+(2.*expect(k2.*du,l1,l2))*N);
        end
        ke = ke + kek2;
    end
    
    if mat.r
        r = evalparampc(mat,'r',PC,elem,xnode,xgauss);
        u = N*qe;
        if ~israndom(u)
            r = expect(r,l1,l2);
            ker = N'*((3.*r.*u.*u)*N);
        else
            ker = N'*(expect(3.*r.*u.*u,l1,l2)*N);
        end
        ke = ke + ker;
    end
end
