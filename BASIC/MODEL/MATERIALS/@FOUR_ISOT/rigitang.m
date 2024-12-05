function ke = rigitang(mat,elem,xnode,xgauss,qe,varargin)
% function ke = rigitang(mat,elem,xnode,xgauss,qe,varargin)

k = evalparam(mat,'k',elem,xnode,xgauss);
B = calc_B(elem,xnode,xgauss);
ke = B'*k*B;

if mat.b || mat.k2 || mat.r || mat.r2 || mat.r3
    N = calc_N(elem,xnode,xgauss);
    
    if mat.b
        b = evalparam(mat,'b',elem,xnode,xgauss);
        ke = ke + (B'*b)*N;
    end
    
    if mat.k2
        k2 = evalparam(mat,'k2',elem,xnode,xgauss);
        u = N*qe;
        du = B*qe;
        ke = ke + B'*k2*((u.*u)*B+(2.*u.*du)*N);
    end
    
    if mat.r
        r = evalparam(mat,'r',elem,xnode,xgauss);
        ke = ke + N'*r*N;
    end
    
    if mat.r2
        r2 = evalparam(mat,'r2',elem,xnode,xgauss);
        u = N*qe;
        ke = ke + N'*r2*(2*(u)*N);
    end
    
    if mat.r3
        r3 = evalparam(mat,'r3',elem,xnode,xgauss);
        u = N*qe;
        ke = ke + N'*r3*(3*(u.*u)*N);
    end
end

switch getdim(elem)
    case 1
        if isparam(mat,'S')
            ke = evalparam(mat,'S',elem,xnode,xgauss); % cross-section area
            ke = ke*S;
        end
    case 2
        if isparam(mat,'DIM3')
            e = evalparam(mat,'DIM3',elem,xnode,xgauss); % thickness
            ke = ke*e;
        end
end
