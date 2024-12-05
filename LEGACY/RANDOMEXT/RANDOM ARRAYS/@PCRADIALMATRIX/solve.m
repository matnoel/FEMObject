function [u,varargout] = solve(A,b,varargin)
% function [u,varargout] = solve(A,b)
% Resolution de Ax=b
% A : PCRADIALMATRIX
% b : PCMATRIX ou PCRADIALMATRIX ou double
%
% function [u,varargout] = solve(A,b,PC)
% PC : POLYCHAOS sur lequel est defini b et la solution x (A peut etre defini sur un autre)
%
% function [u,varargout] = solve(A,b,PC,'norandom',norandom)
% les composantes norandom de x sont deterministes

varargout = cell(1,nargout-1);
callpcmatrix = 1;

if nargout==2
    varargout{1} = 0;
end

if callpcmatrix && isa(A,'PCMATRIX')
    
    if israndom(A)
        A = PCMATRIX(A);
    end
    if israndom(b)
        b = PCMATRIX(b);
    end
    
    [u,varargout{:}] = solve(A,b,varargin{:});
    
else
    
    nosto = getcharin('norandom',varargin);
    
    [A,b,PC] = pcsystemupdate(A,b,varargin{:});
    
    
    n = size(b,1);
    nns = length(nosto);
    sto = setdiff(1:n,nosto);
    ns = length(sto);
    
    if ~israndom(A)
        u = b;
        u.V = A\b.V;
        
        if ~isempty(nosto)
            error('pas programme')
        end
        
    else
        
        if ~isempty(nosto)
            error('non programme')
        end
        
        
        masse = getmasse(A);
        massep = double(masse);
        
        b = double(PCMATRIX(b));
        b1 = b(sto,:);
        B = b1(:);
        
        if iscell(A.V)
            val11 = A.V(sto,sto);
            val11 = getvalue(multimtimes(massep,val11));
            
        else
            val11 = double(A.V(sto,sto));
            val11 = val11*massep';
        end
        
        
        M11 = assemble(val11,ns,ns,length(PC),length(PC));
        M = M11;
        v = zeros(size(M,2),1);
        repzero = find(abs(diag(M))<eps);
        repnonzero = setdiff(1:size(M,2),repzero);
        v(repnonzero) = M(repnonzero,repnonzero)\B(repnonzero);
        v1 = v(1:ns*(length(PC)));
        v2 = v(ns*(length(PC))+1:end);
        u = zeros(n,length(PC));
        u(sto,:) = reshape(v1,ns,length(PC));
        u(nosto,:) = v2*double(one(PC));
        u = PCMATRIX(u,[n,1],PC);
    end
    
    
    
end


function M = assemble(val,n1,n2,p1,p2)
% function M = assemble(val,n1,n2,p1,p2)

if isa(val,'cell')
    val = reshape(val,p1,p2);
    VAL = [];
    I = [];
    J = [];
    for k=1:p1
        for l=1:p2
            [i,j,v] = find(val{k,l});
            i = i+(k-1)*n1;
            j = j+(l-1)*n2;
            I = [I;i];
            J = [J;j];
            VAL = [VAL;v];
        end
    end
    M = sparse(I,J,VAL,n1*p1,n2*p2);
    
else
    [I,J] = find(val);
    [Ji,Jj] = ind2sub([p1,p2],J);
    [Ii,Ij] = ind2sub([n1,n2],I);
    j = (Jj(:)-1)*n2+Ij(:);
    i = (Ji(:)-1)*n1+Ii(:);
    I = sub2ind(size(val),I,J);
    M = sparse(i,j,val(I),n1*p1,n2*p2);
end

return



