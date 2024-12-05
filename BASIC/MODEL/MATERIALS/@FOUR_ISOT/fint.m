function fe = fint(mat,elem,xnode,xgauss,qe,varargin)
% function fe = fint(mat,elem,xnode,xgauss,qe,varargin)

k = evalparam(mat,'k',elem,xnode,xgauss);
B = calc_B(elem,xnode,xgauss);
gradu = B*qe;
fe = B'*(k*gradu);

if mat.b || mat.k2 || mat.r || mat.r2 || mat.r3
    N = calc_N(elem,xnode,xgauss);
    
    if mat.b
        b = evalparam(mat,'b',elem,xnode,xgauss);
        fe = fe + N'*(b'*gradu);
    end
    
    if mat.k2
        k2 = evalparam(mat,'k2',elem,xnode,xgauss);
        u = N*qe;
        fe = fe + B'*(k2*(u.*u)*gradu);
    end
    
    if mat.r
        r = evalparam(mat,'r',elem,xnode,xgauss);
        u = N*qe;
        fe = fe + N'*(r*u);
    end
    
    if mat.r2
        r2 = evalparam(mat,'r2',elem,xnode,xgauss);
        u = N*qe;
        fe = fe + N'*(r2*(u.*u));
    end
    
    if mat.r3
        r3 = evalparam(mat,'r3',elem,xnode,xgauss);
        u = N*qe;
        fe = fe + N'*(r3*(u.*u.*u));
    end
end

switch getdim(elem)
    case 1
        if isparam(mat,'S')
            fe = evalparam(mat,'S',elem,xnode,xgauss); % cross-section area
            fe = fe*S;
        end
    case 2
        if isparam(mat,'DIM3')
            e = evalparam(mat,'DIM3',elem,xnode,xgauss); % thickness
            fe = fe*e;
        end
end
