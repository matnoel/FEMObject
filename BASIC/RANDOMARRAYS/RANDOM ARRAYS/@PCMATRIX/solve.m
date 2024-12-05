function [u,varargout] = solve(A,b,varargin)
% function [u,varargout] = solve(A,b)
% Resolution de Ax=b
% A : PCMATRIX
% b : PCMATRIX ou PCRADIALMATRIX ou double
%
% function [u,varargout] = solve(A,b,PC)
% PC : POLYCHAOS sur lequel est defini b et la solution x (A peut etre defini sur un autre)
%
% function [u,varargout] = solve(A,b,PC,'norandom',norandom)
% les composantes norandom de x sont deterministes
%
% function [u,varargout] = solve(A,b,'singular')
% permet d'eliminer les lignes et colonnes associees a des termes diagonaux
% nuls


nosto = getcharin('norandom',varargin);

[A,b,PC] = pcsystemupdate(A,b,varargin{:});

n = size(b,1);
nns = length(nosto);
sto = setdiff(1:n,nosto);
ns = length(sto);

if nargout==2
    varargout{1} = 0;
end

if ~israndom(A)
    
    b = double(cell2mat(b));
    u = A\b;
    u = PCMATRIX(u,[n,1],PC);
    
    if ~isempty(nosto)
        error('pas programme')
    end
    
else
    
    
    b = double(PCMATRIX(cell2mat(b)));
    
    b1 = b(sto,:);
    
    masse = getmasse(A);
    masseun  =  masse * sparse(double(ones(1,PC))');
    massep = double(masse);
    masseunp = double(masseun);
    
    if iscell(A.MULTIMATRIX)
        val11 = A.MULTIMATRIX(sto,sto);
        val11 = getvalue(multimtimes(massep,val11));
    else
        val11 = double(A.MULTIMATRIX(sto,sto));
        val11 = val11*massep';
    end
    
    M11 = assemble(val11,ns,ns,length(PC),length(PC));
    
    
    if ~isempty(nosto)
        b2 = b(nosto,:)*double(one(PC))';
        B = [b1(:);b2(:)];
        
        if ~iscell(A.MULTIMATRIX)
            val12 = double(A.MULTIMATRIX(sto,nosto));
            val21 = double(A.MULTIMATRIX(nosto,sto));
            val22 = double(A.MULTIMATRIX(nosto,nosto));
            val12 = val12*masseunp';
            val21 = val21*masseunp';
            val22 = val22*sparse(double(one(getPC(A)))');
            
        else
            
            val12 = (A.MULTIMATRIX(sto,nosto));
            val21 = (A.MULTIMATRIX(nosto,sto));
            val22 = (A.MULTIMATRIX(nosto,nosto));
            val12 = getvalue(multimtimes(masseunp,val12));
            val21 = getvalue(multimtimes(masseunp,val21));
            val22 = getvalue(multimtimes(sparse(double(one(getPC(A)))')',val22));
            
        end
        
        M12 = assemble(val12,ns,nns,length(PC),1);
        M21 = assemble(val21,nns,ns,1,length(PC));
        M22 = assemble(val22,nns,nns,1,1);
        
        M = [M11,M12;M21,M22];
        
    else
        B = b1(:);
        M = M11;
    end
    
    v = zeros(size(M,2),1);
    if ischarin('singular',varargin)
        repzero = find(abs(diag(M))<eps);
    else
        repzero = [];
    end
    repnonzero = setdiff(1:size(M,2),repzero);
    
    
    v(repnonzero) = M(repnonzero,repnonzero)\B(repnonzero);
    v1 = v(1:ns*(length(PC)));
    v2 = v(ns*(length(PC))+1:end);
    
    u = zeros(n,length(PC));
    u(sto,:) = reshape(v1,ns,length(PC));
    u(nosto,:) = v2*double(one(PC));
    
    u = PCMATRIX(u,[n,1],PC);
    
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


