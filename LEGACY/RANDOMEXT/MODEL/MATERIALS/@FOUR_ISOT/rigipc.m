function ke = rigipc(mat,elem,xnode,xgauss,PC,varargin)
% function ke = rigipc(mat,elem,xnode,xgauss,PC,varargin)

k = evalparampc(mat,'k',PC,elem,xnode,xgauss);
B = calc_B(elem,xnode,xgauss);
ke = B'*k*B;

if mat.b || mat.k2 || mat.r
    N=calc_N(elem,xnode,xgauss);
    
    if mat.b
        b = evalparampc(mat,'b',PC,elem,xnode,xgauss);
        
        keb = N'*(b'*B);
        
        if isparam(mat,'stabilize') && getparam(mat,'stabilize')>0
            bxi = norm(b);
            he = 2/sum(abs(b'*B/bxi));
            if getparam(mat,'stabilize')==1
                Pe = (bxi*he/2/k);
                tau = he/2./bxi.*(1./tanh(Pe)-1./Pe);
            elseif getparam(mat,'stabilize')==2
                tau=he/2/bxi;
            end
            keb = keb + tau*(B'*b)*(b'*B);
        end
        
        if israndom(ke) && ~israndom(keb)
            keb=one(PC)*keb;
        end
        ke = ke + keb;
    end
    
    if mat.k2 || mat.r
        warning('la matrice de rigidite ne tient pas compte des termes non-lineaires')
    end
    
end

