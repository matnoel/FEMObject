function varargout = eig(A,varargin)
% function e = eig(A)
% function [V,D] = eig(A)
% function [V,D,W] = eig(A)
% function [...] = eig(A,'sort')
% Same function as EIG but for objects of class MYDOUBLEND 
% with an additional option for sorted eigenvalues and eigenvectors

sorted = ischarin('sort',varargin);
varargin = delonlycharin('sort',varargin);

s = size2D(A);
sND = sizeND(A);
pND = prod(sND);

if ischarin('vector',varargin) || (nargout==1 && ~ischarin('matrix',varargin))
    D = zeros([s(1),1,sND]);
else
    D = zeros([s,sND]);
end
if nargout>1
    V = zeros([s,sND]);
end
if nargout>2
    W = zeros([s,sND]);
end
for k=1:pND
    Ak = A.double(:,:,k);
    if nargout==1
        Dk = eig(Ak,varargin{:});
    elseif nargout==2
        [Vk,Dk] = eig(Ak,varargin{:});
    else
        [Vk,Dk,Wk] = eig(Ak,varargin{:});
    end
    if sorted
        if ischarin('vector',varargin) || (nargout==1 && ~ischarin('matrix',varargin))
            [Dk,ind] = sort(Dk);
        else
            [Dk,ind] = sort(diag(Dk));
            Dk = diag(Dk);
        end
    end
    D(:,:,k) = Dk;
    if nargout>1
        if sorted
            V(:,:,k) = Vk(:,ind);
        else
            V(:,:,k) = Vk;
        end
    end
    if nargout>2
        if sorted
            W(:,:,k) = Wk(:,ind);
        else
            W(:,:,k) = Wk;
        end
    end
end

if nargout==1
    varargout{1} = MYDOUBLEND(D);
else
    varargout{1} = MYDOUBLEND(V);
    varargout{2} = MYDOUBLEND(D);
    if nargout>2
         varargout{3} = MYDOUBLEND(W);
    end
end
