function ke = rigiexpect(mat,elem,xnode,xgauss,l1,l2,varargin)
% function ke = rigiexpect(mat,elem,xnode,xgauss,l1,l2,varargin)

k = evalparampc(mat,'k',[],elem,xnode,xgauss);
k = expect(k,l1,l2);
B = calc_B(elem,xnode,xgauss);
ke = B'*k*B;

if mat.b || mat.k2
    N=calc_N(elem,xnode,xgauss);
    
    if mat.b
        b = evalparampc(mat,'b',[],elem,xnode,xgauss);
        b = expect(b,l1,l2);
        ke = ke + (B'*b)*N;
    end
    if mat.k2
        warning('la matrice de rigidite ne tient pas compte des termes non-lineaires')
    end
    
end
