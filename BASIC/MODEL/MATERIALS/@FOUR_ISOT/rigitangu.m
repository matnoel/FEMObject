function keU = rigitangu(mat,elem,xnode,xgauss,qe,Ue,varargin)
% function keU = rigitangu(mat,elem,xnode,xgauss,qe,Ue,varargin)

k = evalparampc(mat,'k',[],elem,xnode,xgauss);
B = calc_B(elem,xnode,xgauss);
BU = B*Ue;
keU = B'*k*BU;

if mat.b || mat.k2 || mat.r
    N = calc_N(elem,xnode,xgauss);
    NU = N*Ue;
    if mat.b
        b = evalparampc(mat,'b',[],elem,xnode,xgauss);
        kebU = (B'*b)*NU;
        keU = keU + kebU;
    end
    
    if mat.k2
        k2 = evalparampc(mat,'k2',[],elem,xnode,xgauss);
        u = N*qe;
        du = B*qe;
        kek2U = B'*((k2.*u.*u)*BU+(2.*k2.*u.*du)*NU);
        keU = keU + kek2U;
    end
    
    if mat.r
        r = evalparampc(mat,'r',[],elem,xnode,xgauss);
        u = N*qe;
        kerU = N'*((3.*r.*u.*u)*NU);
        keU = keU + kerU;
    end
end

switch getdim(elem)
    case 1
        if isparam(mat,'S')
            keU = evalparam(mat,'S',elem,xnode,xgauss); % cross-section area
            keU = keU*S;
        end
    case 2
        if isparam(mat,'DIM3')
            e = evalparam(mat,'DIM3',elem,xnode,xgauss); % thickness
            keU = keU*e;
        end
end
