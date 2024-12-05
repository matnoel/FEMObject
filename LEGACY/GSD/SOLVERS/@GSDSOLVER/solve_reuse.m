function [u,flag,result,bini] = solve_reuse(GSD,A,b,ureuse,varargin)
% function [u,flagreuse,result,bini] = solve_reuse(GSD,A,b,ureuse,varargin)
% resolution de Au=b sur la base reduite issue de reuse
% A et b : PCMATRIX ou PCRADIALMATRIX
% ureuse : pour la reutilisation (MULTIMATRIX ou PCRADIALMATRIX ou double)
% u : solution : PCRADIALMATRIX
% bini = b-A*u (residu)


paramradial.reuse = getparam(GSD,'reuse');
if paramradial.reuse && ~isempty(ureuse)
    if isa(ureuse,'PCRADIALMATRIX')
        if length(ureuse)==0
            paramradial.reuse=false;
        end
    elseif   isa(ureuse,'PCMATRIX')
        fprintf('  on ne peut reutiliser une PCMATRIX -> utiliser PCRADIALMATRIX ou MULTIMATRIX\n')
        paramradial.reuse=false;
        
    elseif ~isa(ureuse,'MULTIMATRIX')
        
        if isa(ureuse,'double') && normest(ureuse)<eps
            fprintf('  vecteur de norme negligeable -> pas d''initialisation\n')
            paramradial.reuse=false;
        end
        
    end
else
    paramradial.reuse = false;
end


if paramradial.reuse
    n=size(A,1);
    
    if isa(ureuse,'PCRADIALMATRIX')
        ureuse = getV(ureuse);
    end
    
    toliter = getparam(GSD,'toliter');
    directsolve = getparam(GSD,'direct');
    
    if directsolve
        localstosolver = @(aU,fU) solve(aU,fU);
    else
        localstosolver = @(aU,fU) cgs(aU,fU,toliter,[],'noupdate');
    end
    
    
    U = double(ureuse);
    fU = U'*b;
    aU = U'*A*U;
    
    
    [l,flagiter] = localstosolver(aU,fU);
    
    result.Rayg = full(expect(fU,l));
    if isa(result.Rayg,'MULTIMATRIX')
        result.Rayg   = double(cell2mat(result.Rayg));
    end
    result.rayg = trace(result.Rayg);
    
    u = PCRADIALMATRIX(U,[n,1],l);
    if nargout>=4
        bini = b - A*u;
    end
    
else
    
    u = zeros(size(A,1),1);
    result.Rayg = 0;
    result.rayg = 0;
    bini = b;
end

flag= ~paramradial.reuse;

